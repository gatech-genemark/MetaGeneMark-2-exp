import os

runningOnPycharm = "PYCHARM_HOSTED" in os.environ
# os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
import matplotlib
# matplotlib.use("pgf")
# print(matplotlib.get_cachedir())
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

if not runningOnPycharm:
    matplotlib.use('Agg')

