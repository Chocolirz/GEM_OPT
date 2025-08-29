import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

params = params = {
    'axes.labelsize': 21,
    'font.size': 20,
    'font.family': 'sans-serif', 
    'font.serif': 'Arial', 
    'legend.fontsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18, 
    'axes.labelpad': 15,
    
    'figure.figsize': [10,8], # value in inches based on dpi of monitor
    'figure.dpi': 105.5, # My monitor has a dpi of around 105.5 px/inch
   
    'axes.grid': False,
    'grid.linestyle': '-',
    'grid.alpha': 0.25,
    'axes.linewidth': 1,
    'figure.constrained_layout.use': True,
   
   
    # Using Paul Tol's notes:
    'axes.prop_cycle': 
        mpl.cycler(color=['#4477aa', # blue
                          '#ee6677', # red/pink
                          '#228833', # green
                          '#aa3377', # purple
                          '#66ccee', # cyan
                          '#ccbb44', # yellow
                          '#bbbbbb', # grey
                          'k']),
       
       # Pick either the cycler above, or the cycler below:
        
       # (mpl.cycler(color=['#4477aa', # blue
       #                     '#ee6677', # red/pink
       #                     '#228833', # green
       #                     '#aa3377', # purple
       #                     '#66ccee', # cyan
       #                     '#ccbb44', # yellow
       #                     '#bbbbbb', # grey
       #                     ]) + 
       #   mpl.cycler(linestyle=['-', # solid
       #                         '--', # dashed
       #                         ':', # dotted
       #                         '-.', # dash dot
       #                         (0, (3, 1, 1, 1, 1, 1)), # narrow dash dot dot
       #                         (0, (1, 2, 5, 2, 5, 2)), # dash dash dot
       #                         (0, (5, 2.5, 1, 2.5, 1, 2.5)), # dash dot dot
       #                         ])), 
       
    'lines.linewidth': 2.5,
   
    'image.cmap': 'turbo',#'PuBuGn',
   
    'legend.frameon': False,
   
    'legend.edgecolor': 'white',
   
    'legend.framealpha': 0.5,
} 

plt.rcParams.update(params)

log_tick_format = ticker.FuncFormatter(lambda y,pos: f'{y:.{int(np.maximum(-np.log10(y),0)):1d}f}')

ratio_aligned = np.array([[0.0, 0.0, 0.36357315014349223], 
                 [0.0, 0.5, 0.72313741230891], 
                 [0.0, 1.0, 1.0720667456041368], 
                 [0.0, 1.5, 1.401782129179304], 
                 [0.0, 2.0, 1.7005637798829885], 
                 [0.0, 2.5, 1.9499812586714873], 
                 [0.0, 3.0, 2.1124738400787457], 
                 [0.0, 3.5, 2.232164468719458], 
                 [0.0, 4.0, 2.325139280933423], 
                 [0.5, 0.0, 2.687162107776381], 
                 [0.5, 0.5, 3.045121251651384], 
                 [0.5, 1.0, 3.392185219293731], 
                 [0.5, 1.5, 3.720503590182157], 
                 [0.5, 2.0, 4.0165420306556054], 
                 [0.5, 2.5, 4.261578073213777], 
                 [0.5, 3.0, 4.422292557199408], 
                 [0.5, 3.5, 4.541133084831325], 
                 [0.5, 4.0, 4.633817800833074], 
                 [1.0, 0.0, 4.987317843789462], 
                 [1.0, 0.5, 5.336641412617146], 
                 [1.0, 1.0, 5.674977980344922], 
                 [1.0, 1.5, 5.991996268864941], 
                 [1.0, 2.0, 6.276609117603689], 
                 [1.0, 2.5, 6.502933622820072], 
                 [1.0, 3.0, 6.6536994996035], 
                 [1.0, 3.5, 6.767510222077461], 
                 [1.0, 4.0, 6.856785997561392], [1.5, 0.0, 7.1928584651467675], [1.5, 0.5, 7.52562166141697], 
                 [1.5, 1.0, 7.844569800641119], [1.5, 1.5, 8.141100675036784], [1.5, 2.0, 8.398263108410763], 
                 [1.5, 2.5, 8.581598896463511], [1.5, 3.0, 8.717498286863481], [1.5, 3.5, 8.82376151859181], 
                 [1.5, 4.0, 8.908238954490828], [2.0, 0.0, 9.216210676699957], [2.0, 0.5, 9.5198102201481], 
                 [2.0, 1.0, 9.809151472760993], [2.0, 1.5, 10.068786140720514], [2.0, 2.0, 10.26820506580835], 
                 [2.0, 2.5, 10.419974757966925], [2.0, 3.0, 10.540739022335678], [2.0, 3.5, 10.637569613790204], 
                 [2.0, 4.0, 10.716244692852884], [2.5, 0.0, 10.978894123476703], [2.5, 0.5, 11.235477979436478], 
                 [2.5, 1.0, 11.469781702239814], [2.5, 1.5, 11.658338818359885], [2.5, 2.0, 11.812808403867223], 
                 [2.5, 2.5, 11.939987488861837], [2.5, 3.0, 12.045315889223334], [2.5, 3.5, 12.132466311247136], 
                 [2.5, 4.0, 12.204644175322857], [3.0, 0.0, 12.383245878126194], [3.0, 0.5, 12.55666576293386], 
                 [3.0, 1.0, 12.71843107062806], [3.0, 1.5, 12.862491078870956], [3.0, 2.0, 12.98750339692888], 
                 [3.0, 2.5, 13.094128679342978], [3.0, 3.0, 13.186309050722079], [3.0, 3.5, 13.263819245978834], 
                 [3.0, 4.0, 13.329694301923798], [3.5, 0.0, 13.464668176105214], [3.5, 0.5, 13.597737751351469], 
                 [3.5, 1.0, 13.722785484276534], [3.5, 1.5, 13.836799409644547], [3.5, 2.0, 13.940038765436578], 
                 [3.5, 2.5, 14.031035161281698], [3.5, 3.0, 14.109959096116844], [3.5, 3.5, 14.178105548589702], 
                 [3.5, 4.0, 14.236877443478031], [4.0, 0.0, 14.34399810144012], [4.0, 0.5, 14.449325636617566], 
                 [4.0, 1.0, 14.550001971177663], [4.0, 1.5, 14.643586132518568], [4.0, 2.0, 14.73085531547307], 
                 [4.0, 2.5, 14.809319520343116], [4.0, 3.0, 14.877602152785563], [4.0, 3.5, 14.937375094460073], 
                 [4.0, 4.0, 14.98950540398284]])

ratio_misaligned = np.array([[0.0, 0.0, 0.36291026619473415], [0.0, 0.5, 0.7218141929184028], 
                             [0.0, 1.0, 1.0701045173206087], [0.0, 1.5, 1.399131818802771], 
                             [0.0, 2.0, 1.6972082388873626], [0.0, 2.5, 1.9458344635670952], 
                             [0.0, 3.0, 2.10853200331235], [0.0, 3.5, 2.2283864084179235], [0.0, 4.0, 2.3214505492272086], 
                             [0.5, 0.0, 2.6828373238235494], [0.5, 0.5, 3.040093187845833], [0.5, 1.0, 3.3864604115632932], 
                             [0.5, 1.5, 3.7141015527881813], [0.5, 2.0, 4.009371355609591], [0.5, 2.5, 4.253588311760447], 
                             [0.5, 3.0, 4.414524592212236], [0.5, 3.5, 4.533466957854146], [0.5, 4.0, 4.626251759151834], 
                             [1.0, 0.0, 4.9790132478986315], [1.0, 0.5, 5.327702055483825], [1.0, 1.0, 5.665322078225918], 
                             [1.0, 1.5, 5.981498145492813], [1.0, 2.0, 6.265391387682247], [1.0, 2.5, 6.491023114738329], 
                             [1.0, 3.0, 6.641940154367703], [1.0, 3.5, 6.755863849067516], [1.0, 4.0, 6.845221949041498], 
                             [1.5, 0.0, 7.180637452996422], [1.5, 0.5, 7.512571259578695], [1.5, 1.0, 7.830687344950116], 
                             [1.5, 1.5, 8.126574835678735], [1.5, 2.0, 8.38271245247966], [1.5, 2.5, 8.566200144277298], 
                             [1.5, 3.0, 8.70227650731612], [1.5, 3.5, 8.808655783396055], [1.5, 4.0, 8.893184133238595], 
                             [2.0, 0.0, 9.200254875652346], [2.0, 0.5, 9.50318116778235], [2.0, 1.0, 9.791734466257715], 
                             [2.0, 1.5, 10.05015169093941], [2.0, 2.0, 10.249697276341607], [2.0, 2.5, 10.40163083349298], 
                             [2.0, 3.0, 10.52246942248548], [2.0, 3.5, 10.619383315590252], [2.0, 4.0, 10.698134626401592], 
                             [2.5, 0.0, 10.959784954766977], [2.5, 0.5, 11.215182824188284], [2.5, 1.0, 11.448708692770765],
                             [2.5, 1.5, 11.637598699864533], [2.5, 2.0, 11.792050960580333], [2.5, 2.5, 11.919270419239067], 
                             [2.5, 3.0, 12.024712279725668], [2.5, 3.5, 12.111949957232994], [2.5, 4.0, 12.184200013326077], 
                             [3.0, 0.0, 12.362995488418612], [3.0, 0.5, 12.536674772016116], [3.0, 1.0, 12.6984543252184], 
                             [3.0, 1.5, 12.84238270945007], [3.0, 2.0, 12.967579121848495], [3.0, 2.5, 13.07432047690264], 
                             [3.0, 3.0, 13.166574111636264], [3.0, 3.5, 13.244155146943092], [3.0, 4.0, 13.310109584463135], 
                             [3.5, 0.0, 13.445022716110799], [3.5, 0.5, 13.578024756709604], [3.5, 1.0, 13.703148764860915], 
                             [3.5, 1.5, 13.81734726658955], [3.5, 2.0, 13.920601265884908], [3.5, 2.5, 14.01166915116009], 
                             [3.5, 3.0, 14.09069220885265], [3.5, 3.5, 14.158914619918697], [3.5, 4.0, 14.217749421949666], 
                             [4.0, 0.0, 14.324938915412092], [4.0, 0.5, 14.430345901763719], [4.0, 1.0, 14.53107429772788], 
                             [4.0, 1.5, 14.624719260326525], [4.0, 2.0, 14.712114564304844], [4.0, 2.5, 14.790664569368523], 
                             [4.0, 3.0, 14.859013687867023], [4.0, 3.5, 14.91885317100001], [4.0, 4.0, 14.971037435074349]])

ratio_new_a = ratio_aligned[:,2]
for i in range(1,len(ratio_new_a)):
    ratio_new_a[-i] = ratio_new_a[-i]-ratio_new_a[-i-1]

ratio_new = ratio_misaligned[:,2]
for i in range(1,len(ratio_new)):
    ratio_new[-i] = ratio_new[-i]-ratio_new[-i-1]

PlotScatter = False
if PlotScatter:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')


    ratio = np.array(ratio_aligned)  # make it a numpy array for convenience
    ax.scatter(ratio[:,0], ratio[:,1], ratio_new_a, marker='o')

    ratio = np.array(ratio_misaligned)  # make it a numpy array for convenience
    ax.scatter(ratio[:,0], ratio[:,1], ratio_new, marker='o')

    ax.set_xlabel(r'$x$ [cm]')
    ax.set_ylabel(r'$y$ [cm]')
    ax.set_zlabel(r'Collection ratio')

    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(ratio[:,0], ratio[:,1], ratio_new_a - ratio_new, marker='o')
    ax.set_xlabel(r'$x$ [cm]')
    ax.set_ylabel(r'$y$ [cm]')
    ax.set_zlabel(r'Collection ratio diff.')

    plt.show()

PlotSurface = False
if PlotSurface:
    # Extract x, y, z
    x = ratio_aligned[:,0]
    y = ratio_aligned[:,1]
    z = ratio_new_a # or ratio_new

    # Find unique values (since it's on a grid)
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    # Reshape z into a 2D grid
    X, Y = np.meshgrid(y_unique, x_unique)   # careful: order
    Z = z.reshape(len(x_unique), len(y_unique))

    # Plot surface
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='plasma', edgecolor='none')

    ax.set_xlabel(r'$x$ [cm]')
    ax.set_ylabel(r'$y$ [cm]')
    ax.set_zlabel(r'Collection ratio diff.')
    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.show()

PlotDiff = False
if PlotDiff:
    # Extract x, y, z
    x = ratio_aligned[:,0]
    y = ratio_aligned[:,1]
    z = ratio_new_a - ratio_new

    # Unique grid values
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    # Reshape z into grid
    X, Y = np.meshgrid(y_unique, x_unique)  # note: order!
    Z = z.reshape(len(x_unique), len(y_unique))

    # Plot 2D heatmap
    plt.figure(figsize=(7,6))
    c = plt.pcolormesh(X, Y, Z, cmap="plasma", shading='auto')  # shading='auto' prevents artifacts
    plt.colorbar(c, label=r'Collection ratio diff.')

    plt.xlabel(r'$y$ [cm]')
    plt.ylabel(r'$x$ [cm]')

    plt.show()

PlotR = True
if PlotR:
    # Extract x, y, z
    x = ratio_aligned[:,0]
    y = ratio_aligned[:,1]
    r = np.sqrt(x**2 + y**2)
    z = ratio_new_a - ratio_new

    plt.figure(figsize=(8,8))
    plt.xlabel(r'$R$ [cm]')
    plt.ylabel(r'Collection ratio diff.')
    plt.plot(r, z, 'o')
    plt.grid(True)
    plt.show()