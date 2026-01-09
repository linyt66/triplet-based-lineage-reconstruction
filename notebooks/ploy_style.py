
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
# -------------------------- 1. 顶刊风格全局配置（通用，无需修改）--------------------------
# 建议放在所有绘图代码最前面
import matplotlib.pyplot as plt

plt.rcParams.update({
    # —— 输出与字体嵌入（矢量优先，可在 Illustrator/AI 中编辑）——
    'savefig.format': 'pdf',          # Nature 推荐提交 PDF（矢量）
    'figure.dpi': 300,                # 预览分辨率
    'savefig.dpi': 600,               # 栅格元素（若有）更清晰
    'pdf.fonttype': 42,               # TrueType，便于后期编辑
    'ps.fonttype': 42,
    'svg.fonttype': 'none',           # SVG 保留文字，不转曲线

    # —— 版式与尺寸
    # （单栏单图（3，2），双栏双图(5.5, 4)，双栏三图(7.5, 1.5)，可根据宽度高度稍微调整）——
    'figure.figsize': (3, 2),     # 默认单栏尺寸
    'figure.facecolor': 'white',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'savefig.transparent': False,     # Journal 通常要求白底
    'figure.constrained_layout.use': False,

    # —— 字体（Nature：无衬线，常用 Arial/Helvetica）——
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],  # 加 DejaVu 兜底
    'axes.unicode_minus': True,       # 使用真正的“−”号
    'text.usetex': False,             # 非必须不启用 LaTeX

    # —— 数学字体与指数记法（与无衬线主体保持一致）——
    'mathtext.fontset': 'dejavusans',
    'axes.formatter.use_mathtext': True,
    'axes.formatter.limits': (-3, 3), # 仅在 <1e-3 或 >1e3 时切换科学计数

    # —— 字号（印刷 6.5–7 pt 常用；一致、紧凑）——
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5,
    'legend.title_fontsize': 7,

    # —— 坐标轴与刻度（四边框 + 内向刻度 + 启用次刻度）——
    'axes.spines.top': True,  # 设置为false就是去掉图框中上面的线
    'axes.spines.right': True, # 设置为false就是去掉图框中右面的线
    'axes.linewidth': 0.6,            # 细边框（~0.5–0.7 pt）
    'axes.labelpad': 2.0,

    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,

    'xtick.minor.visible': True,  # 显示 x 轴的次刻度线
    'ytick.minor.visible': True,  # 显示 y 轴的次刻度线

    'xtick.major.size': 4,
    'xtick.major.width': 0.6,
    'xtick.minor.size': 2,
    'xtick.minor.width': 0.4,
    'xtick.major.pad': 2,

    'ytick.major.size': 4,
    'ytick.major.width': 0.6,
    'ytick.minor.size': 2,
    'ytick.minor.width': 0.4,
    'ytick.major.pad': 2,

    # —— 线型、误差棒与图例 —— 
    'lines.linewidth': 0.8,           # 细线，避免过粗
    'lines.markersize': 3.0,
    'errorbar.capsize': 2.0,

    'legend.frameon': False,
    'legend.handlelength': 2.0,
    'legend.handleheight': 0.7,
    'legend.borderaxespad': 1,
    'legend.columnspacing': 2,
})

def log_format(x, pos):
    if x > 0:
        return f'{int(np.log10(x))}'  # 将10^x的标签简化为x
    else:
        return f'{x}'

# 顶刊配色（通用）
colors = {
    'blue': '#1f77b4', 'orange': '#ff7f0e', 'green': '#2ca02c',
    'red': '#d62728', 'purple': '#9467bd', 'dark_gray': '#595959'
}
