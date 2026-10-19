"""Rebuild final five-method figures, standalone HTML and printable PDF from frozen tables."""
import base64,csv,html,json,math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.lib.pagesizes import A4

ROOT=Path(__file__).resolve().parents[1]; PAPER=ROOT/'paper'; FIG=PAPER/'figures'
METHODS=['MSIpath','MSIsensor-pro','MSIsensor2','mSINGS','MANTIS2']
DRAW_ORDER=METHODS[1:]+['MSIpath']
COLORS=dict(zip(METHODS,['#174c83','#32936f','#9467bd','#de8f05','#d65f5f']))
def line(ax,x,y,method):
 return ax.plot(x,y,label=method,color=COLORS[method],marker='D' if method=='MSIpath' else 'o',
    markerfacecolor='white' if method=='MSIpath' else COLORS[method],
    linewidth=2.3 if method=='MSIpath' else 1.3,markersize=6 if method=='MSIpath' else 4,
    zorder=10 if method=='MSIpath' else 3)
def rows(path):
 with path.open(encoding='utf8') as f:return list(csv.DictReader(f,delimiter='\t'))
DATA={name:rows(PAPER/'data'/name/'performance_metrics.tsv') for name in ['end_to_end','mechanism','robustness']}
RES=rows(PAPER/'data/resources/resource_comparison.tsv'); MEAS=rows(PAPER/'data/resources/measurements.tsv')
def select(name,method,scope='panel200',**kw):
 found=[r for r in DATA[name] if r['method']==method and r['scope']==scope and all(str(r[k])==str(v) for k,v in kw.items())]
 assert len(found)<=1,(name,method,kw)
 return found[0] if found else None
def pct(x):return f'{100*float(x):g}%'
def rate(r):return pct(r['matched_rate']) if r else '—'
def fnum(v,n=3):return f'{float(v):.{n}f}'
def finish(fig,name):
 fig.tight_layout();fig.savefig(FIG/(name+'.png'),dpi=220,bbox_inches='tight')
 fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def figures():
 FIG.mkdir(parents=True,exist_ok=True)
 for scope in ['panel200','model_common']:
  fig,axes=plt.subplots(2,2,figsize=(10.4,7.2));conditions=['B10_P20','B10_P60','B30_P20','B30_P60']
  for ax,c in zip(axes.flat,conditions):
   rr=[select('end_to_end',m,scope,condition=c) for m in METHODS]; rr=[r for r in rr if r]
   x=np.arange(len(rr));v=np.array([float(r['matched_rate']) for r in rr]);lo=np.array([float(r['matched_CI_low']) for r in rr]);hi=np.array([float(r['matched_CI_high']) for r in rr])
   ax.bar(x,v,color=[COLORS[r['method']] for r in rr],width=.6)
   ax.errorbar(x,v,yerr=[v-lo,hi-v],fmt='none',ecolor='#303030',capsize=3)
   ax.set_xticks(x,[r['method'] for r in rr],rotation=20,ha='right');ax.set_ylim(0,1.07);ax.set_ylabel('Sensitivity')
   b,p=c.split('_');ax.set_title(f'b={int(b[1:])}%, p={int(p[1:])}% | 30x | n=200');ax.grid(axis='y',alpha=.15)
  finish(fig,'01_end_to_end_'+scope)
  from sklearn.metrics import roc_curve,precision_recall_curve,roc_auc_score,average_precision_score
  scored=rows(PAPER/'data/end_to_end/sample_scores.tsv')
  for kind in ['ROC','PR']:
   fig,axes=plt.subplots(2,2,figsize=(10.4,7.2))
   for ax,c in zip(axes.flat,conditions):
    for m in DRAW_ORDER:
     rr=[r for r in scored if r['role']=='test' and r['scope']==scope and r['method']==m and r['condition'] in ['N',c] and r['score']!='']
     if not rr:continue
     y=np.array([int(r['label']) for r in rr]);v=np.array([float(r['score']) for r in rr])
     if len(set(y))<2:continue
     if kind=='ROC':
      x,z,_=roc_curve(y,v);metric=roc_auc_score(y,v)
     else:
      z,x,_=precision_recall_curve(y,v);metric=average_precision_score(y,v)
     ax.plot(x,z,color=COLORS[m],linewidth=2.3 if m=='MSIpath' else 1.2,zorder=10 if m=='MSIpath' else 3,label=f'{m}: {metric:.3f} (n={len(rr)})')
    bpart,ppart=c.split('_');ax.set_title(f'b={int(bpart[1:])}%, p={int(ppart[1:])}% | 30x')
    ax.set_xlabel('False-positive rate' if kind=='ROC' else 'Recall');ax.set_ylabel('True-positive rate' if kind=='ROC' else 'Precision')
    ax.set_xlim(-.02,1.02);ax.set_ylim(-.02,1.02);ax.legend(fontsize=7,loc='lower right' if kind=='ROC' else 'lower left');ax.grid(alpha=.15)
   finish(fig,'01_'+kind.lower()+'_'+scope)
  fig,axes=plt.subplots(1,2,figsize=(12,4.6))
  negatives=['N','allele_imbalance','rare_in_panel','novel_normal']
  for m in DRAW_ORDER:
   rr=[select('mechanism',m,scope,condition=c) for c in negatives]
   if any(r is None for r in rr):continue
   line(axes[0],range(4),[float(r['matched_rate']) for r in rr],m)
   line(axes[1],range(4),[float(r['matched_unresolved']) for r in rr],m)
  for ax,title in zip(axes,['False-positive rate','Unresolved rate']):
   ax.set_xticks(range(4),['Common\nnormal','Allele\nimbalance','Rare in\npanel','Novel\nnormal']);ax.set_ylim(-.03,1.07);ax.set_title(title+' | n=50 per condition');ax.grid(alpha=.15)
  axes[1].legend(fontsize=9);finish(fig,'02_mechanism_'+scope)
  fig,axes=plt.subplots(2,2,figsize=(10.4,7.2))
  for ax,depth in zip(axes.flat,[10,20,30,60]):
   for m in DRAW_ORDER:
    rr=[r for r in DATA['robustness'] if r['scope']==scope and r['method']==m and int(r['depth'])==depth and r['direction']=='mixed' and r['endpoint']=='sensitivity']
    rr.sort(key=lambda r:float(r['p']))
    if rr:line(ax,[100*float(r['p']) for r in rr],[float(r['matched_rate']) for r in rr],m)
   ax.set_title(f'{depth}x | b=30% | n=50');ax.set_xlabel('Tumor purity p (%)');ax.set_ylabel('Sensitivity');ax.set_ylim(-.03,1.07);ax.grid(alpha=.15)
  axes[0,0].legend(fontsize=8);finish(fig,'03_robustness_'+scope)
  fig,axes=plt.subplots(1,2,figsize=(11,4.3))
  for ax,p in zip(axes,[.2,.8]):
   for m in DRAW_ORDER:
    rr=[]
    for direction in ['mixed','expansion','contraction']:
     found=[r for r in DATA['robustness'] if r['scope']==scope and r['method']==m and int(r['depth'])==30 and float(r['p'])==p and r['direction']==direction]
     if len(found)==1:rr.append(found[0])
    if len(rr)==3:line(ax,range(3),[float(r['matched_rate']) for r in rr],m)
   ax.set_xticks(range(3),['Mixed','Expansion','Contraction']);ax.set_title(f'30x | p={int(p*100)}% | b=30%');ax.set_ylabel('Sensitivity');ax.set_ylim(-.03,1.07);ax.grid(alpha=.15)
  axes[1].legend(fontsize=8);finish(fig,'04_direction_'+scope)
  rr=[r for r in RES if r['scope']==scope]
  fig,axes=plt.subplots(1,3,figsize=(13,4.3))
  for ax,field,title in zip(axes,['wall_seconds','peak_tree_PSS_MiB','peak_tree_RSS_MiB'],['BAM-to-call time (s)','Peak tree PSS (MiB)','Peak tree RSS (MiB)']):
   samples=[[float(v[field]) for v in MEAS if v['scope']==scope and v['method']==r['method']] for r in rr]
   ax.boxplot(samples,showfliers=True);ax.set_xticks(range(1,len(rr)+1),[r['method'] for r in rr]);ax.tick_params(axis='x',rotation=30);ax.set_title(title);ax.grid(axis='y',alpha=.15)
   if field=='wall_seconds':ax.set_yscale('log')
  fig.suptitle('Same 10 samples | 3 repetitions | 1 CPU');finish(fig,'05_resources_'+scope)

class Document:
 def __init__(self):
  self.story=[];self.html=[];self.text=[]
  styles=getSampleStyleSheet()
  font='CJK'
  candidates=[Path('C:/Windows/Fonts/simhei.ttf'),Path('/mnt/c/Windows/Fonts/simhei.ttf'),Path('/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc')]
  for p in candidates:
   if p.exists():pdfmetrics.registerFont(TTFont(font,str(p)));break
  else:pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'));font='STSong-Light'
  self.font=font
  self.body=ParagraphStyle('BodyCN',fontName=font,fontSize=9.2,leading=14,spaceAfter=7,wordWrap='CJK',textColor=colors.HexColor('#213547'))
  self.heading=ParagraphStyle('HeadCN',parent=self.body,fontSize=17,leading=24,spaceAfter=13,textColor=colors.HexColor('#174c83'))
  self.small=ParagraphStyle('SmallCN',parent=self.body,fontSize=7.6,leading=11,spaceAfter=4)
 def head(self,title):
  self.story.append(Paragraph(html.escape(title),self.heading));self.html.append('<h2>'+html.escape(title)+'</h2>');self.text.append('\n'+title+'\n')
 def para(self,text):
  self.story.append(Paragraph(html.escape(text),self.body));self.html.append('<p>'+html.escape(text)+'</p>');self.text.append(text)
 def tab(self,headers,rs,widths=None):
  data=[[Paragraph(html.escape(str(x)),self.small) for x in row] for row in [headers]+rs]
  widths=widths or [515/len(headers)]*len(headers)
  t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
  t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9f1f8')),('VALIGN',(0,0),(-1,-1),'TOP'),
   ('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#174c83')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f7f9fb')]),
   ('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
  self.story.extend([t,Spacer(1,10)])
  self.html.append('<table><thead><tr>'+''.join('<th>'+html.escape(str(c))+'</th>' for c in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(c))+'</td>' for c in r)+'</tr>' for r in rs)+'</tbody></table>')
  self.text.append('\t'.join(headers));self.text.extend('\t'.join(map(str,r)) for r in rs)
 def fig(self,name,width=515):
  p=FIG/(name+'.png')
  from PIL import Image as PILImage
  with PILImage.open(p) as im:w,h=im.size
  self.story.append(Image(str(p),width=width,height=width*h/w));self.story.append(Spacer(1,8))
  self.html.append('<img src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'">')
 def page(self):self.story.append(PageBreak());self.html.append('<hr>')
 def save(self):
  out=PAPER/'output/pdf';out.mkdir(parents=True,exist_ok=True)
  def footer(canvas,doc):
   canvas.setFont(self.font,8);canvas.setFillColor(colors.HexColor('#63758a'))
   canvas.drawString(40,23,'MSIpath | RECOMB 实验汇总 | 2026-10-13');canvas.drawRightString(555,23,str(doc.page))
  SimpleDocTemplate(str(out/'MSIpath_final_experiment_report.pdf'),pagesize=A4,leftMargin=40,rightMargin=40,topMargin=38,bottomMargin=42,
    title='MSIpath - Final experimental report',author='MSIpath research project').build(self.story,onFirstPage=footer,onLaterPages=footer)
  css='body{max-width:1100px;margin:40px auto;padding:0 24px;font:16px/1.7 system-ui;color:#213547}h2{color:#174c83}table{border-collapse:collapse;width:100%;font-size:14px;margin:20px 0}td,th{padding:9px;border-bottom:1px solid #dce4ec;text-align:left}th{background:#e9f1f8}img{width:100%;height:auto}hr{margin:40px 0;border:0;border-top:1px solid #dce4ec}@media print{hr{break-after:page}body{margin:0}img,table{break-inside:avoid}}'
  (PAPER/'MSIpath_final_experiment_report.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>MSIpath 最终实验报告</title><style>'+css+'</style><body>'+''.join(self.html)+'</body></html>',encoding='utf8')
  (PAPER/'SUMMARY_ZH.txt').write_text('\n\n'.join(self.text)+'\n',encoding='utf8')

def main():
 figures();d=Document()
 d.head('MSIpath：RECOMB 投稿实验汇总')
 d.para('版本：MSIpath 1.0.0。模型、结果命名及论文图表已统一；比较对象为 MSIsensor-pro、MSIsensor2、mSINGS、MANTIS2。四项模拟实验及五软件统一资源测量均已完成。')
 d.para('主要结果：MSIpath在30倍深度、30%异常位点、20%肿瘤纯度时检出率为83.5%，普通正常误报率为2%；同深度纯扩张场景检出率82%，高于MSIsensor-pro的6%。但10%纯度检出较弱，罕见及库外正常压力测试误报严重，运行时间和内存也高于对比工具。')
 d.tab(['实验','供体 / 测试实例','设计'],[
  ['端到端','200 / 1,000','30x；b=10%、30%；p=20%、60%；另有普通正常'],
  ['机制','50 / 250','普通正常、罕见正常、库外正常、库外体细胞异常、80:20等位基因不平衡'],
  ['鲁棒性','50 / 1,200','10/20/30/60x；b=30%；p=0/10/20/40/80%；扩张/收缩对照'],
  ['资源','5 / 10，每软件3次','15共同位点五软件；200位点四软件；270次新进程测量']],widths=[76,119,320])
 d.para('模拟采用chr19和chr22目标区域、200个STR位点、PE150及ART PCR-free测序误差。p表示肿瘤纯度，b表示异常位点比例；模型假设二倍体、拷贝数中性、克隆性单等位基因改变，改变的分子比例为p/2。不同实验复用供体，不把实例总数当成独立患者数。')
 d.para('主判读采用同深度独立正常队列校准，显著性阈值0.05；结果表另保留默认阈值及固定30x校准。各软件实际误报率不同。15共同位点与200位点结果分开报告；MANTIS2使用额外配对正常信息。')
 d.para('本报告用于整理投稿实验证据，不将模拟结果表述为临床性能，也不将开发阶段复用数据称为独立盲测。')
 d.page();d.head('1. 端到端检测性能')
 conditions=['B10_P20','B10_P60','B30_P20','B30_P60']
 for scope,title in [('panel200','200位点面板'),('model_common','15共同位点：五软件比较')]:
  d.para(title+'。下表为检出率，末列为普通正常误报率；每个条件200名供体。')
  rs=[]
  for m in METHODS:
   rs.append([m]+[rate(select('end_to_end',m,scope,condition=c)) for c in conditions]+[rate(select('end_to_end',m,scope,condition='N'))])
  d.tab(['软件','b10/p20','b10/p60','b30/p20','b30/p60','正常误报'],rs,widths=[105,82,82,82,82,82])
 d.para('MSIsensor2现有官方模型只覆盖15个共同位点，200位点表中的空缺不代表0检出率。MSIpath在200位点四个阳性场景的AUROC依次为0.7931、0.9932、0.9878、1.0000。83.5%检出率的95% Wilson区间为77.7%–88.0%；2%普通正常误报率的区间为0.8%–5.0%。')
 d.para('未判定样本保留在检出/误报的分母中；较低误报不必然意味着更多正确阴性。各条件的可评分率、未判定率和正确MSS率见随附完整表。MSIsensor2和mSINGS在当前小位点集上的表现不能直接外推到其原始训练面板。')
 d.page();d.head('1. 端到端：条件间变化与置信区间')
 d.fig('01_end_to_end_panel200')
 d.para('误差线为95% Wilson区间。MSIpath在b=10%、p=60%条件优于MSIsensor-pro（89.5%对77%），但b=10%、p=20%仍低于后者（14%对49%）；不能概括为所有条件均领先。')
 d.para('AUROC/AUPRC的阴性对照为同深度普通正常样本。此类指标并未纳入下一节的罕见/库外正常压力测试，不能用高AUROC替代对困难正常背景的检验。')
 d.page();d.head('2. 机制实验：遗传背景与体细胞异常')
 cs=['N','allele_imbalance','rare_in_panel','novel_normal','novel_somatic'];rs=[]
 for m in METHODS:
  if m=='MSIsensor2':continue
  rs.append([m]+[rate(select('mechanism',m,condition=c)) for c in cs])
 d.tab(['软件','普通正常\n误报','不平衡\n误报','罕见正常\n误报','库外正常\n误报','体细胞异常\n检出'],rs,widths=[105,82,82,82,82,82])
 d.fig('02_mechanism_panel200')
 d.para('每种条件50个实例。MSIpath在普通正常及不平衡条件下误报为2%和0%，但库内罕见正常、库外正常均为100%；库外体细胞异常检出100%。这些结果说明异常识别与异常来源区分是不同问题，当前模型尚未可靠解决后者。')
 d.para('罕见/库外正常是保留SNP背景、定向编辑约40个STR的合成压力测试，不能将100%外推为真实一般人群误报率。MSIsensor2的机制结果单独保留于15共同位点表和图。')
 d.page();d.head('3. 深度与肿瘤纯度鲁棒性')
 d.fig('03_robustness_panel200')
 d.para('固定b=30%，每个条件50名供体。MSIpath在30x、p=10/20/40/80%时检出率为12/84/100/100%；60x时为22/100/100/100%。10%纯度仍明显弱于MSIsensor-pro（30x为82%，60x为90%）。')
 d.para('MSIpath在p=20%时，10/20/30/60x检出率为16/34/84/100%；相应普通正常误报率为2/0/2/2%，均可判读。MSIsensor-pro在10x条件全部未判定，图中0检出不能解释为有效阴性判读。')
 d.page();d.head('3. 长度变化方向对照')
 d.fig('04_direction_panel200')
 d.tab(['30x条件','MSIpath','MSIsensor-pro','解释'],[
  ['p=20%，纯扩张','82%','6%','扩张识别优势'],['p=20%，纯收缩','82%','100%','收缩仍有差距'],
  ['p=80%，纯扩张','100%','4%','高纯度扩张也可检出'],['p=80%，纯收缩','100%','100%','均达到高检出']],widths=[145,85,105,180])
 d.para('纯扩张与纯收缩对照使用匹配的目标位点、单倍型及改变长度绝对值。混合方向组并非所有情况下都与两种强制方向组具有完全相同的异常集合，图中连线用于比较条件，不把它解释为连续变量趋势。')
 d.page();d.head('4. 五软件统一资源基准')
 d.para('同一批10个样本来自5名供体，正常与p=20%、b=30%的MSI样本各5个；每软件3次重复。固定1 CPU、1线程、串行运行。机器为Intel Core Ultra 9 285，WSL2环境。下表均为30次执行的中位数；内存列是每次执行峰值的中位数。')
 for scope,title in [('model_common','15共同位点：五软件'),('panel200','200位点：四软件')]:
  d.para(title)
  rr=[r for r in RES if r['scope']==scope]
  d.tab(['软件','耗时(s)','CPU(s)','PSS(MiB)','RSS(MiB)','可判读样本'],[
   [r['method'],fnum(r['wall_seconds_median']),fnum(r['cpu_seconds_median']),fnum(r['peak_tree_PSS_MiB_median'],1),fnum(r['peak_tree_RSS_MiB_median'],1),r['scorable_samples']+'/10'] for r in rr],widths=[115,80,80,80,80,80])
 d.para('时间覆盖进程启动、载入预建模型、BAM证据提取、推断、最终校准判读及写出。mSINGS包含名称适配、索引、mpileup、过滤、analyzer和count；MANTIS2包含配对正常处理。排除模拟、公共比对、一次性参考/基线训练和输出归档。')
 d.para('MSIpath的耗时与内存均高于对比工具，不能声称计算效率领先。15位点范围下mSINGS仅5/10、MANTIS2仅9/10可判读；失败于证据不足的执行仍保留资源记录，不能把更快结束当成准确率优势。')
 d.page();d.head('4. 资源分布与计量边界')
 d.fig('05_resources_model_common')
 d.para('图中每软件30次执行，时间为对数轴。重复运行不是30名独立供体，不据此做患者层面的显著性检验。OS缓存未清空、执行顺序固定种子随机化。')
 d.para('PSS按共享比例计入内存，RSS可能重复计入共享页。20ms采样可能漏掉短时峰值；原始表同时给出GNU time最大单进程RSS高水位。共同轻量Python封装也计入资源。全部270次执行的得分与可用位点数均与冻结实验一致。')
 d.para('既有分阶段资源测试显示200位点HMM/图构建约41–45秒，而单等位基因似然约0.40–0.75秒；进一步提速应优先针对路径推断。分阶段测试使用32/64/128候选压力配置，与正式Top20配置的统一资源测试不能直接相加。发布版重构保留算法，但封装和文件写出成本尚未单独重测。')
 d.page();d.head('5. 消融诊断、结论及投稿边界')
 d.para('MSIpath内部参考诊断：罕见正常干预位点的名义显著比例，个体化混合为64.45%，仅MAP基因型为60.10%，全panel人群模型为37.20%。库外正常分别为76.75%、65.45%、80.10%；库外体细胞异常分别为63.40%、43.25%、71.00%。这些是位点层面的正常秩诊断，不是样本检出率或BH校正结果；全panel同时改变参考支持与权重，不能当作纯先验消融。')
 d.para('可支持的结论：个体化路径背景下的单等位基因模型能整合多位点证据；在本模拟的中等纯度和纯扩张条件有较强检出能力；普通正常与不平衡条件下误报较低。不能支持的结论：所有场景优于对照、解决10%纯度检出、解决所有罕见生殖系混淆、计算效率领先或获得临床验证。')
 d.para('独立性：模型开发的小测试供体包含在后续复用队列中，其他场景也经过探索查看。冻结参数后扩展测试不等于新的盲测。正式泛化结论仍需供体独立的验证集；40名正常校准样本的最小秩p值为1/41，阈值稳定性仍需更大正常队列验证。')
 d.para('可复现交付：代码包仅有MSIpath默认模型入口，附带研究参考面板、冻结正常误差参数、校准库、两例证据演示、单元测试及数值一致性验证。所有原始实验保留；paper/data中仅转换方法显示名及纯度单位，分数、判读、区间和失败场景未改动。')
 d.para('表格与图：端到端、机制、鲁棒性各有200位点和15共同位点版本；资源有五软件共同位点与四软件200位点版本。原始CSV/TSV、图像PNG及矢量PDF均随包附带。')
 d.para('参考：Huang et al. ART: a next-generation sequencing read simulator. Bioinformatics (2012), doi:10.1093/bioinformatics/btr708；EnsembleTR群体STR资源，https://pmc.ncbi.nlm.nih.gov/articles/PMC10593948/。软件来源与版本见docs/SOURCES.md及validation/environment.json。')
 d.page();d.head('附图：ROC曲线（200位点）');d.fig('01_roc_panel200')
 d.para('图例为AUROC及可评分实例数。每个条件对比相同深度的普通正常；曲线只使用有分数的实例，样本级检出率仍使用完整分母。15共同位点版本随图像文件附带，其可评分实例数可能因软件而异。')
 d.page();d.head('附图：精确率-召回率曲线（200位点）');d.fig('01_pr_panel200')
 d.para('图例为平均精确率（AUPRC）及可评分实例数。模拟的正负样本构成影响AUPRC，不代表临床筛查患病率下的阳性预测值。所有主比较阈值由正常校准决定，没有使用测试ROC选择部署阈值。')
 d.save()
 print(PAPER/'output/pdf/MSIpath_final_experiment_report.pdf')
if __name__=='__main__':main()
