import pygmo as pg
from pylab import plt, np
def my_plot_non_dominated_fronts(points, marker='o', comp=[0, 1], up_to_rank_no=None):
    # We plot
    fronts, _, _, _ = pg.fast_non_dominated_sorting(points)

    # We define the colors of the fronts (grayscale from black to white)
    if up_to_rank_no is None:
        cl = list(zip(np.linspace(0.1, 0.9, len(fronts)),
                      np.linspace(0.1, 0.9, len(fronts)),
                      np.linspace(0.1, 0.9, len(fronts))))
    else:
        cl = list(zip(np.linspace(0.1, 0.9, up_to_rank_no),
                      np.linspace(0.1, 0.9, up_to_rank_no),
                      np.linspace(0.1, 0.9, up_to_rank_no)))

    fig, ax = plt.subplots()

    count = 0
    for ndr, front in enumerate(fronts):
        count += 1
        # We plot the points
        for idx in front:
            ax.plot(points[idx][comp[0]], points[idx][
                comp[1]], marker=marker, color=cl[ndr])
        # We plot the fronts
        # Frist compute the points coordinates
        x = [points[idx][comp[0]] for idx in front]
        y = [points[idx][comp[1]] for idx in front]
        # Then sort them by the first objective
        tmp = [(a, b) for a, b in zip(x, y)]
        tmp = sorted(tmp, key=lambda k: k[0])
        # Now plot using step
        ax.step([c[0] for c in tmp], [c[1]
                                      for c in tmp], color=cl[ndr], where='post')
        if up_to_rank_no is None:
            pass
        else:
            if count >= up_to_rank_no:
                break

    return ax
def my_2p5d_plot_non_dominated_fronts(points, marker='o', comp=[0, 1], up_to_rank_no=1, no_text=True):
    # from pylab import mpl
    # mpl.rcParams['font.family'] = ['Times New Roman']
    # mpl.rcParams['font.size'] = 16.0

    full_comp = [0, 1, 2]
    full_comp.remove(comp[0])
    full_comp.remove(comp[1])
    z_comp = full_comp[0]

    # We plot
    # fronts, dl, dc, ndr = pg.fast_non_dominated_sorting(points)
    fronts, _, _, _= pg.fast_non_dominated_sorting(points)

    # We define the colors of the fronts (grayscale from black to white)
    cl = list(zip(np.linspace(0.9, 0.1, len(fronts)),
                  np.linspace(0.9, 0.1, len(fronts)),
                  np.linspace(0.9, 0.1, len(fronts))))

    fig, ax = plt.subplots(constrained_layout=False)
    plt.subplots_adjust(left=None, bottom=None, right=0.85, top=None, wspace=None, hspace=None)

    count = 0
    for ndr, front in enumerate(fronts):
        count += 1

        # Frist compute the points coordinates
        x_scale = 1
        y_scale = 1
        z_scale = 1
        if comp[0] == 1: # efficency
            x_scale = 100
        if comp[1] == 1: # efficency
            y_scale = 100
        if z_comp == 1: # efficency
            z_scale = 100
        x = [points[idx][comp[0]]*x_scale for idx in front]
        y = [points[idx][comp[1]]*y_scale for idx in front]
        z = [points[idx][z_comp] *z_scale for idx in front]

        # # We plot the points
        # for idx in front:
        #     ax.plot(points[idx][comp[0]], points[idx][comp[1]], marker=marker, color=cl[ndr])

        # Then sort them by the first objective
        tmp = [(a, b, c) for a, b, c in zip(x, y, z)]
        tmp = sorted(tmp, key=lambda k: k[0])
        # Now plot using step
        ax.step([coords[0] for coords in tmp], 
                [coords[1] for coords in tmp], color=cl[ndr], where='post')

        # Now add color according to the value of the z-axis variable usign scatter
        scatter_handle = ax.scatter(x, y, c=z, alpha=0.5, cmap='Spectral', marker=marker, zorder=99) #'viridis'
        # color bar
        cbar_ax = fig.add_axes([0.875, 0.15, 0.02, 0.7])
        cbar_ax.get_yaxis().labelpad = 10
        clb = fig.colorbar(scatter_handle, cax=cbar_ax)
        if z_comp == 0:
            z_label = r'$-\rm {TRV}$ [Nm/m^3]'
            z_text  = '%.0f'
        elif z_comp == 1:
            z_label = r'$-\eta$ [%]'
            z_text  = '%.1f'
        elif z_comp == 2:
            z_label = r'$O_C$ [1]'
            z_text  = '%.1f'
        clb.ax.set_ylabel(z_label, rotation=270)


        if z_comp == 2: # when OC as z-axis
            print('-----------------------------------------------------')
            print('-----------------------------------------------------')
            print('-----------------------------------------------------')
            # Add index next to the points
            for x_coord, y_coord, z_coord, idx in zip(x, y, z, front):
                if no_text:
                    pass
                else:
                    ax.annotate( z_text%(z_coord) + ' #%d'%(idx), (x_coord, y_coord) )
        else:
            # text next scatter showing the value of the 3rd objective
            for i, val in enumerate(z):
                if no_text:
                    pass
                else:
                    ax.annotate( z_text%(val), (x[i], y[i]) )

        # refine the plotting
        if comp[0] == 0:
            ax.set_xlabel(r'$-\rm {TRV}$ [Nm/m^3]')
        elif comp[0] == 1:
            ax.set_xlabel(r'$-\eta$ [%]')
        elif comp[0] == 2:
            ax.set_xlabel(r'$O_C$ [1]')

        if comp[1] == 0:
            ax.set_ylabel(r'$-\rm {TRV}$ [Nm/m^3]')
        elif comp[1] == 1:
            ax.set_ylabel(r'$-\eta$ [%]')
        elif comp[1] == 2:
            ax.set_ylabel(r'$O_C$ [1]')
        ax.grid()

        # plot up to which domination rank?
        if up_to_rank_no is None:
            pass
        else:
            if count >= up_to_rank_no:
                break
    # Y730
    # fig.savefig(r'C:\Users\horyc\Desktop/'+ '2p5D-%d%d.png'%(comp[0],comp[1]), dpi=300)
    return ax
def my_3d_plot_non_dominated_fronts(pop, paretoPoints, fea_config_dict, az=180, comp=[0, 1, 2], plot_option=1):
    """
    Plots solutions to the DTLZ problems in three dimensions. The Pareto Front is also
    visualized if the problem id is 2,3 or 4.
    Args:
        pop (:class:`~pygmo.population`): population of solutions to a dtlz problem
        az (``float``): angle of view on which the 3d-plot is created
        comp (``list``): indexes the fitness dimension for x,y and z axis in that order
    Returns:
        ``matplotlib.axes.Axes``: the current ``matplotlib.axes.Axes`` instance on the current figure
    Raises:
        ValueError: if *pop* does not contain a DTLZ problem (veryfied by its name only) or if *comp* is not of length 3
    Examples:
        >>> import pygmo as pg
        >>> udp = pg.dtlz(prob_id = 1, fdim =3, dim = 5)
        >>> pop = pg.population(udp, 40)
        >>> udp.plot(pop) # doctest: +SKIP
    """
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    import numpy as np
    # from pylab import mpl
    # mpl.rcParams['font.family'] = ['Times New Roman']
    # mpl.rcParams['font.size'] = 16.0

    # if (pop.problem.get_name()[:-1] != "DTLZ"):
    #     raise(ValueError, "The problem seems not to be from the DTLZ suite")

    if (len(comp) != 3):
        raise(ValueError, "The kwarg *comp* needs to contain exactly 3 elements (ids for the x,y and z axis)")

    # Create a new figure
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')

    # plot the points
    fits = np.transpose(pop.get_f())
    try:
        pass
        # ax.plot(fits[comp[0]], fits[comp[1]], fits[comp[2]], 'ro')
    except IndexError:
        print('Error. Please choose correct fitness dimensions for printing!')

    if False:
        # Plot pareto front for dtlz 1
        if plot_option==1: # (pop.problem.get_name()[-1] in ["1"]):

            X, Y = np.meshgrid(np.linspace(0, 0.5, 100), np.linspace(0, 0.5, 100))
            Z = - X - Y + 0.5
            # remove points not in the simplex
            for i in range(100):
                for j in range(100):
                    if X[i, j] < 0 or Y[i, j] < 0 or Z[i, j] < 0:
                        Z[i, j] = float('nan')

            ax.set_xlim(0, 1.)
            ax.set_ylim(0, 1.)
            ax.set_zlim(0, 1.)

            ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
            plt.plot([0, 0.5], [0.5, 0], [0, 0])

        # Plot pareto fronts for dtlz 2,3,4
        if plot_option == 2: # (pop.problem.get_name()[-1] in ["2", "3", "4"]):
            # plot the wireframe of the known optimal pareto front
            thetas = np.linspace(0, (np.pi / 2.0), 30)
            # gammas = np.linspace(-np.pi / 4, np.pi / 4, 30)
            gammas = np.linspace(0, (np.pi / 2.0), 30)

            x_frame = np.outer(np.cos(thetas), np.cos(gammas))
            y_frame = np.outer(np.cos(thetas), np.sin(gammas))
            z_frame = np.outer(np.sin(thetas), np.ones(np.size(gammas)))

            ax.set_autoscalex_on(False)
            ax.set_autoscaley_on(False)
            ax.set_autoscalez_on(False)

            ax.set_xlim(0, 1.8)
            ax.set_ylim(0, 1.8)
            ax.set_zlim(0, 1.8)

            ax.plot_wireframe(x_frame, y_frame, z_frame)

    # https://stackoverflow.com/questions/37000488/how-to-plot-multi-objectives-pareto-frontier-with-deap-in-python
    # def simple_cull(inputPoints, dominates):
    #     paretoPoints = set()
    #     candidateRowNr = 0
    #     dominatedPoints = set()
    #     while True:
    #         candidateRow = inputPoints[candidateRowNr]
    #         inputPoints.remove(candidateRow)
    #         rowNr = 0
    #         nonDominated = True
    #         while len(inputPoints) != 0 and rowNr < len(inputPoints):
    #             row = inputPoints[rowNr]
    #             if dominates(candidateRow, row):
    #                 # If it is worse on all features remove the row from the array
    #                 inputPoints.remove(row)
    #                 dominatedPoints.add(tuple(row))
    #             elif dominates(row, candidateRow):
    #                 nonDominated = False
    #                 dominatedPoints.add(tuple(candidateRow))
    #                 rowNr += 1
    #             else:
    #                 rowNr += 1

    #         if nonDominated:
    #             # add the non-dominated point to the Pareto frontier
    #             paretoPoints.add(tuple(candidateRow))

    #         if len(inputPoints) == 0:
    #             break
    #     return paretoPoints, dominatedPoints
    # def dominates(row, candidateRow):
    #     return sum([row[x] >= candidateRow[x] for x in range(len(row))]) == len(row)  
    # import random
    # print(inputPoints)
    # inputPoints = [[random.randint(70,100) for i in range(3)] for j in range(500)]
    # print(inputPoints)
    # quit()
    # inputPoints = [(x,y,z) for x,y,z in zip(fits[comp[0]], fits[comp[1]], fits[comp[2]])]
    # paretoPoints, dominatedPoints = simple_cull(inputPoints, dominates)
    x = [coords[0]/1000 for coords in paretoPoints]
    y = [coords[1] for coords in paretoPoints]
    z = [coords[2] for coords in paretoPoints]

    # from surface_fitting import surface_fitting
    # surface_fitting(x,y,z)
    # quit()

    if False:
        pass
    else:
        import pandas as pd
        from pylab import cm
        print(dir(cm))
        df = pd.DataFrame({'x': x, 'y': y, 'z': z})

        # 只有plot_trisurf这一个函数，输入是三个以为序列的，其他都要meshgrid得到二维数组的(即ndim=2的数组) 
        # # https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html
            # surf = ax.plot_trisurf(df.x, df.y, df.z, cmap=cm.magma, linewidth=0.1, edgecolor='none')
            # surf = ax.plot_trisurf(x, y, z, cmap='viridis', edgecolor='none')
        # surf = ax.plot_trisurf(df.x, df.y, df.z, cmap=cm.magma, linewidth=0.1)
        surf = ax.plot_trisurf(df.x, df.y*100, df.z, cmap=cm.Spectral, linewidth=0.1)        

        with open('./%s_PF_points.txt'%(fea_config_dict['run_folder'][:-1]), 'w') as f:
            f.write('TRV,eta,OC\n')
            f.writelines(['%g,%g,%g\n'%(a,b,c) for a,b,c in zip(df.x, df.y*100, df.z)])
        # quit()

        fig.colorbar(surf, shrink=0.5, aspect=5)

        ax.set_xlabel(' \n$\\rm -TRV$ [$\\rm kNm/m^3$]')
        ax.set_ylabel(' \n$-\\eta$ [%]')
        ax.set_yticks(np.arange(-96,-93.5,0.5))
        ax.set_zlabel(r'$O_C$ [1]')

        # Try to export data from plot_trisurf # https://github.com/WoLpH/numpy-stl/issues/19
        # print(surf.get_vector())

        # plt.savefig('./plots/avgErrs_vs_C_andgamma_type_%s.png'%(k))
        # plt.show()

        # # rotate the axes and update
        # for angle in range(0, 360):
        #     ax.view_init(30, angle)
        #     plt.draw()
        #     plt.pause(.001)

    ax.view_init(azim=245, elev=15) 
    # fig.tight_layout()

    # Y730
    # fig.savefig(r'C:\Users\horyc\Desktop/3D-plot.png', dpi=300, layout='tight')

    # ax.view_init(azim=az)
    # ax.set_xlim(0, 1.)
    # ax.set_ylim(0, 1.)
    # ax.set_zlim(0, 10.)
    return ax
def my_plot(fits, vectors, ndf):
    plt.rcParams['mathtext.fontset'] = 'stix' # 'cm'
    plt.rcParams["font.family"] = "Times New Roman"
    if True:
        # for fit in fits:
        #     print(fit)
        # ax = pg.plot_non_dominated_fronts(fits)

        # ax = my_plot_non_dominated_fronts(fits, comp=[0,1], marker='o', up_to_rank_no=3)        
        # ax = my_plot_non_dominated_fronts(fits, comp=[0,2], marker='o', up_to_rank_no=3)
        # ax = my_plot_non_dominated_fronts(fits, comp=[1,2], marker='o', up_to_rank_no=3)

        pass

        ax = my_2p5d_plot_non_dominated_fronts(fits, comp=[0,1], marker='o', up_to_rank_no=1)
        # # for studying LSA population (whether or not the optimal is on Rank 1 Pareto Front)
        # x = fits[0][0]/1e3
        # y = fits[0][1]
        # z = fits[0][2]
        # ax.plot(x, y, color='k', marker='s')
        # ax.annotate(r'$x_{\rm optm}$', xy=(x, y), xytext=(x+1, y+0.0005), arrowprops=dict(facecolor='black', shrink=0.05),)
        ax = my_2p5d_plot_non_dominated_fronts(fits, comp=[0,2], marker='o', up_to_rank_no=1)
        ax = my_2p5d_plot_non_dominated_fronts(fits, comp=[1,2], marker='o', up_to_rank_no=1)

    else:
        # Obselete. Use ax.step instead to plot
        print('Valid for 2D objective function space only for now.')
        fig, axes = plt.subplots(ncols=2, nrows=1, dpi=150, facecolor='w', edgecolor='k');
        ax = axes[0]

        for index, xy in enumerate(vectors):
            ax.plot(*xy, 's', color='0')
            ax.text(*xy, '#%d'%(index), color='0')
        ax.title.set_text("Decision Vector Space")

        ax = axes[1]

        for index, front in enumerate(ndf):
            print('Rank/Tier', index, front)
            the_color = '%g'%(index/len(ndf))
            for individual in front: # individual is integer here
                ax.plot(*fits[individual], 'o',                 color=the_color)
                ax.text(*fits[individual], '#%d'%(individual),  color=the_color)

            front = front.tolist()
            front.sort(key = lambda individual: fits[individual][1]) # sort by y-axis value # individual is integer here
            for individual_A, individual_B in zip(front[:-1], front[1:]): # front is already sorted
                fits_A = fits[individual_A]
                fits_B = fits[individual_B]
                ax.plot([fits_A[0], fits_B[0]], 
                        [fits_B[1], fits_B[1]], color=the_color, lw=0.25) # 取高的点的Y
                ax.plot([fits_A[0], fits_A[0]], 
                        [fits_A[1], fits_B[1]], color=the_color, lw=0.25) # 取低的点的X
        ax.title.set_text("Pareto Front | Objective Function Space")
def my_print(ad, pop, _):
    # ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(fits)
    # extract and print non-dominated fronts
    # - ndf (list of 1D NumPy int array): the non dominated fronts
    # - dl (list of 1D NumPy int array): the domination list
    # - dc (1D NumPy int array): the domination count
    # - ndr (1D NumPy int array): the non domination ranks
    fits, vectors = pop.get_f(), pop.get_x()
    ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(fits)

    with open(ad.solver.output_dir+'MOO_log.txt', 'a', encoding='utf-8') as fname:
        print('-'*40, 'Generation:', _, file=fname)
        for rank_minus_1, front in enumerate(ndf):
            print('Rank/Tier', rank_minus_1+1, front, file=fname)
        index = 0
        for domination_list, domination_count, non_domination_rank in zip(dl, dc, ndr):
            print('Individual #%d\t'%(index), 'Belong to Rank #%d\t'%(non_domination_rank), 'Dominating', domination_count, 'and they are', domination_list, file=fname)
            index += 1

        # print(fits, vectors, ndf)
        print(pop, file=fname)


def learn_about_the_archive(prob, swarm_data, popsize, fea_config_dict, len_s01=None, len_s02=None, bool_plot_and_show=False):
    number_of_chromosome = len(swarm_data)
    print('Archive size:', number_of_chromosome)
    # for el in swarm_data:
    #     print('\t', el)

    pop_archive = pg.population(prob, size=number_of_chromosome)
    for i in range(number_of_chromosome):
        pop_archive.set_xf(i, swarm_data[i][:-3], swarm_data[i][-3:])

    sorted_index = pg.sort_population_mo(points=pop_archive.get_f())
    print('Sorted by domination rank and crowding distance:', len(sorted_index))
    print('\t', sorted_index)

    # 这段代码对于重建种群来说不是必须的，单单sort_population_mo（包含fast_non_dominated_sorting和crowding_distance）就够了，
    # 只是我想看看，具体的crowding_distance是多少，然后我想知道排在前面的多少个是属于domination rank 1的。
    if True:
        fits, vectors = pop_archive.get_f(), pop_archive.get_x()
        ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(fits)

        ind1, ind2 = 0, 0
        for rank_minus_1, front in enumerate(ndf):

            ind2 += len(front)
            sorted_index_at_this_front = sorted_index[ind1:ind2]
            fits_at_this_front = [fits[point] for point in sorted_index_at_this_front]

            # Rank 1 Pareto Front
            if ind1 == 0:
                rank1_ParetoPoints = fits_at_this_front
                if len(front) < popsize:
                    print('There are not enough chromosomes (%d) belonging to domination rank 1 (the best Pareto front).\nWill use rank 2 or lower to reach popsize of %d.'%(len(front), popsize))

            # this crwdsit should be already sorted as well
            if len(fits_at_this_front) >= 2: # or else error:  A non dominated front must contain at least two points: 1 detected.
                crwdst = pg.crowding_distance(fits_at_this_front)
            else:
                print('A non dominated front must contain at least two points: 1 detected.')
                crwdst = [999999]


            print('\nRank/Tier', rank_minus_1+1, 'chromosome count:', len(front), len(sorted_index_at_this_front))
            # print('\t', sorted_index_at_this_front.tolist())
            # print('\t', crwdst)
            print('\tindex in pop\t|\tcrowding distance')
            for index, cd in zip(sorted_index_at_this_front, crwdst):
                print('\t', index, '\t\t\t|', cd)
            ind1 = ind2

    sorted_vectors = [vectors[index].tolist() for index in sorted_index]
    sorted_fits    = [fits[index].tolist() for index in sorted_index]

    if bool_plot_and_show:
        my_plot(pop_archive.get_f(), pop_archive.get_x(), ndf)
        my_3d_plot_non_dominated_fronts(pop_archive, rank1_ParetoPoints, fea_config_dict, plot_option=1)
        plt.show()

    swarm_data_on_pareto_front = [design_parameters_denorm + fits for design_parameters_denorm, fits in zip(sorted_vectors, sorted_fits)]
    return swarm_data_on_pareto_front

def pyx_draw_model(im):
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "Times New Roman"

    myfontsize = 13.5
    plt.rcParams.update({'font.size': myfontsize})

    # # 示意图而已，改改尺寸吧
    # im.Radius_OuterStatorYoke -= 37
    # im.Radius_InnerStatorYoke -= 20
    # im.Radius_Shaft += 20
    # # im.Location_RotorBarCenter2 += 5 # this will change the shape of rotor slot

    import VanGogh
    vg = VanGogh.VanGogh_pyPlotter(im, VanGogh.CUSTOM)
    vg.draw_model()

    # PyX
    import pyx
    vg.tikz.c = pyx.canvas.canvas() # clear the canvas because we want to redraw 90 deg with the data vg.tikz.track_path
    from copy import deepcopy
    def pyx_draw_path(vg, path, sign=1):
        if len(path) == 4:
            vg.tikz.draw_line(path[:2], path[2:4], untrack=True)
        else:
            vg.tikz.draw_arc(path[:2], path[2:4], path[4:6], relangle=sign*path[6], untrack=True)
    def rotate(_, x, y):
        return np.cos(_)*x + np.sin(_)*y, -np.sin(_)*x + np.cos(_)*y
    def is_at_stator(im, path):
        return np.sqrt(path[0]**2 + path[1]**2) > im.Radius_OuterRotor + 0.5*im.Length_AirGap
    
    for path in (vg.tikz.track_path): # track_path is passed by reference and is changed by mirror
        path_mirror = deepcopy(path)
        # for mirror copy (along x-axis)
        path_mirror[1] = path[1]*-1
        path_mirror[3] = path[3]*-1

        # rotate path and plot
        if is_at_stator(im, path):
            Q = im.Qs
        else:
            Q = im.Qr
        _ = 2*np.pi/Q
        path[0], path[1] = rotate(0.5*np.pi - 0.5*_, path[0], path[1])
        path[2], path[3] = rotate(0.5*np.pi - 0.5*_, path[2], path[3])
        pyx_draw_path(vg, path, sign=1)

        path_mirror[0], path_mirror[1] = rotate(0.5*np.pi - 0.5*_, path_mirror[0], path_mirror[1])
        path_mirror[2], path_mirror[3] = rotate(0.5*np.pi - 0.5*_, path_mirror[2], path_mirror[3])
        pyx_draw_path(vg, path_mirror, sign=-1)

        # 注意，所有 tack_path 中的 path 都已经转动了90度了！
        # for mirror copy (along y-axis)
        path[0] *= -1
        path[2] *= -1
        pyx_draw_path(vg, path, sign=-1)

        path_mirror[0] *= -1
        path_mirror[2] *= -1
        pyx_draw_path(vg, path_mirror, sign=1)


    # # 整体转动90度。
    # for path in vg.tikz.track_path:
    #     if is_at_stator(im, path):
    #         Q = im.Qs
    #     else:
    #         Q = im.Qr
    #     _ = 2*np.pi/Q
    #     path[0], path[1] = rotate(0.5*np.pi - _, path[0], path[1])
    #     path[2], path[3] = rotate(0.5*np.pi - _, path[2], path[3])
    #     pyx_draw_path(vg, path)
    # track_path_backup = deepcopy(vg.tikz.track_path)

    # # Rotate Copy
    # for path in deepcopy(vg.tikz.track_path):
    #     if is_at_stator(im, path):
    #         Q = im.Qs
    #     else:
    #         Q = im.Qr
    #     _ = 2*np.pi/Q
    #     path[0], path[1] = rotate(_, path[0], path[1])
    #     path[2], path[3] = rotate(_, path[2], path[3])
    #     pyx_draw_path(vg, path)

    # # Rotate Copy
    # for path in (vg.tikz.track_path):
    #     # if np.sqrt(path[0]**2 + path[1]**2) > im.Radius_OuterRotor + 0.5*im.Length_AirGap:
    #     if is_at_stator(im, path):
    #         Q = im.Qs
    #     else:
    #         Q = im.Qr
    #     _ = 2*np.pi/Q
    #     path[0], path[1] = rotate(_, path[0], path[1])
    #     path[2], path[3] = rotate(_, path[2], path[3])
    #     pyx_draw_path(vg, path, sign=-1)

    vg.tikz.c.writePDFfile("selected_otimal_design%s"%(im.ID))
    # vg.tikz.c.writeEPSfile("pyx_output")
    print('Write to pdf file: selected_otimal_design%s.pdf.'%(im.ID))
    quit()

