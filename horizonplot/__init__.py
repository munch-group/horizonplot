import time
import warnings
from itertools import chain
# import gc

import numpy as np
import pandas as pd

# # Make inline plots vector graphics instead of raster graphics
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('pdf', 'svg')

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from math import isclose, floor, log10

def horizon(row, i, cut):
    """
    Compute the values for the three 
    positive and negative intervals.
    """
    val = getattr(row, i)

    if np.isnan(val):
        for i in range(8):
            yield 0

        # for nan color
        yield cut
    else:
        if val < 0:
            for i in range(4):
                yield 0

        val = abs(val)
        for i in range(3):
            yield min(cut, val)
            val = max(0, val-cut)
        yield int(not isclose(val, 0, abs_tol=1e-8)) * cut

        if val >= 0:
            for i in range(4):
                yield 0

        # for nan color
        yield 0

def chrom_sort(item):
    """
    Sorts in a meaningful way for chromosomes.
    """
    if item.startswith('chr'):
        item = item[3:]
    if item.isdigit():
        return item.zfill(3)
    else:
        return item

def round_to_1_signif(x):
    """
    Rounds to first significant digit.
    """
    return round(x, -int(floor(log10(abs(x)))))
    
def horizonplot(df, key, width, cut='fixed', start='start', col='chrom', row='pop', pop_sorting=None, size=0.5, aspect=40):
    """
    Horizon bar plot made allowing multiple chromosomes and multiple samples.
    """

    pop, chrom = row, col
        
    # set cut if not set
    if cut is 'fixed':
        cut = np.max([np.max(df[key]), np.max(-df[key])]) / 3

    # FIXME: cut should be either:
    # value applied to all
    # 'fixed' (default) computing a shared cutoff for all rows
    # 'adaptive' computing a cutoff that fits each row.

    # make the data frame to plot
    row_iter = df.itertuples()
    col_iterators = zip(*(horizon(row, key, cut) for row in row_iter))
    col_names = ['yp1', 'yp2', 'yp3', 'yp4', 
                 'yn1', 'yn2', 'yn3', 'yn4', 'nan']

    df2 = (df[[key, start, chrom, pop]]
           .assign(**dict(zip(col_names, col_iterators)))
          )

    df3 = pd.DataFrame(dict((col, list(chain.from_iterable(zip(df2[col].values, df2[col].values)))) for col in df2))
    #df3[start] = list(df3[start].values[1:]) + [df3[start].values[-1] + width]
    df3[start] = df3.groupby(row).start.apply(lambda sr: pd.concat([sr.iloc[1:], pd.Series([df3[start].iloc[-1] + width])])).values

    df2 = df3

    # # instead group df2 by col, row and add bounds to each group
    # def add_boundries(df):
    #     df = df.copy()
    #     begin, end = df.iloc[0], df.iloc[-1]
    #     begin[col_names], end[col_names] = 0, 0
    #     df = pd.concat([begin.to_frame().T, df, end.to_frame().T], ignore_index=True)
    #     print(df.head(1))
    #     print(df.tail(1))
    #     print()
    #     return df
    # df2 = df2.groupby([col, row]).apply(add_boundries)

    # chromosome names
    chrom_names = list(df.groupby(chrom).groups.keys())
    sorted_chrom_names = sorted(chrom_names, key=chrom_sort)
    
    if pop_sorting is None:
        pop_sorting = sorted(set(df.reset_index()[pop]))
    
    # number of populations
    nr_pop = len(df.groupby(pop).groups)
    
    # sizes of chromosomes
    chrom_sizes = list()
    for chrom_name in sorted_chrom_names:
        chrom_subset = df.loc[df[chrom] == chrom_name]
        est_chrom_len = np.max(chrom_subset.start) + width
        chrom_sizes.append(est_chrom_len)
        
    # relative width of each plot facet 
    # (using lengths of chromosomes)
    facet_widths_ratios = chrom_sizes# * nr_pop

    # make the plot
    with sns.axes_style("ticks"):

        # ingore UserWarning from seaborn that tight_layout is not applied
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # make the facet grid
            g = sns.FacetGrid(df2, 
                            col=chrom, 
                            row=pop,
                            # sharex=False,
                            sharex=True,
                            sharey=True,
                            # margin_titles=True,
                            height=size, 
                            aspect=aspect,
                            col_order=sorted_chrom_names,
                            row_order=pop_sorting,                      
                            gridspec_kws={'hspace':0.0, 
                                            "width_ratios": facet_widths_ratios}
                            )

            # plot colors
            colours = sns.color_palette("Blues", 3) + ['midnightblue'] + \
                    sns.color_palette("Reds", 3) + ['darkred'] + ['lightgrey']

            # first y tick
            ytic1 = round_to_1_signif(cut / 3)

            for col_name, colour in zip(col_names, colours):
                plt.setp(g.fig.texts, text="") # hack to make y facet labels align...
                # map barplots to each facet
                g.map(plt.fill_between, 
                    start, 
                    col_name, 
                    y2=0,
                    color=colour,
                    linewidth=0,
                    capstyle='butt')

            # g.set_ylabels('')

            def add_pop_labels(pop_label, **kwargs):
                # only rightmosts facets:
                ax = plt.gca()
                if ax.get_position().x1 == max(x.get_position().x1 for x in g.axes.flat):
                    p = pop_label.reset_index(drop=True)[0]
                    plt.annotate(p, xy=(1.005 , 0.5), xycoords='axes fraction', ha='left', size=8)

            g.map(add_pop_labels, pop)

            def add_chrom_labels(chrom_label, **kwargs):
                # only topmost facets:
                ax = plt.gca()
                if ax.get_position().y1 == max(x.get_position().y1 for x in g.axes.flat):
                    p = chrom_label.reset_index(drop=True)[0]
                    plt.annotate(p, xy=(0.5 , 1.1), xycoords='axes fraction', ha='center', size=8)

            g.map(add_chrom_labels, chrom)

            for arr in g.axes:
                for ax, max_val in zip(arr, facet_widths_ratios):
                    ax.set_xlim(0, max_val+1)
                    ax.set_ylim(0, cut)
                    ax.set(xlabel='', ylabel='')
                    ax.set(xticks=np.arange(0, max_val, round_to_1_signif(max_val) / 10))
    #                ax.set_yticks([ytic1, ytic1*2, ytic1*3])
                    ax.set_yticks([ytic1*1.5])
                    g.set_titles('', '')
                
            # remove top and right frame
            sns.despine()

            plt.subplots_adjust(right=0.95)
            
            return g.fig

if __name__ == "__main__":

    n = 1500
    pops = 50

    df = pd.DataFrame({'chrom': ['chr1']*pops*n,
                    'pop': [x for y in ([chr(65+i)]*n for i in range(pops)) for x in y],
                    'start': list(range(1*n)) * pops, 
#                    'pi': list(np.sin(np.linspace(-np.pi, 10*np.pi, 1*n))+0.1) * pops })
                     'pi': np.add(list(np.sin(np.linspace(-np.pi, 10*np.pi, 1*n))) * pops, np.random.random(n*pops)) })
 
    df.loc[(df.start > 40) & (df.start < 60), 'pi'] = np.nan

    g = sns.FacetGrid(df, 
                        col='chrom', 
                        row='pop',
                        sharex=True)
    g.map(plt.plot, 'start', 'pi')

    plt.savefig('tmp2.pdf')

    fig = horizonplot(df, 'pi', width=1, col='chrom', row='pop', size=0.3, aspect=100)
    plt.savefig('tmp.pdf')
    # plt.close(fig)  # close to allow garbage collection, also suppresses inline plot
    # #         gc.collect()