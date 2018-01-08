import time
# import gc

import numpy as np
import pandas as pd

# # Make inline plots vector graphics instead of raster graphics
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('pdf', 'svg')

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
        
class Timer(object):
    
    def __init__(self, verbose=True):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        if self.verbose:
            print('elapsed time: {} secs'.format(self.secs))


def horizon_plot(df, key, width, cut='fixed', start='start', chrom='chrom', pop='pop', pop_sorting=None):
    """
    Horizon bar plot made allowing multiple chromosomes and multiple samples.
    """
    
    from math import isclose, floor, log10
    
    def horizon(row, i, cut):
        """
        Compute the values for the three 
        positive and negative intervals.
        """
        val = getattr(row, i)

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

    def chrom_sort(item):
        """
        Sorts in a meaningful way for chromosomes.
        """
        if item.startswith('chr'):
            item = item[3:]
        if item.isdigit():
            return int(item)
        else:
            return item

    def round_to_1_signif(x):
        """
        Rounds to first significant digit.
        """
        return round(x, -int(floor(log10(abs(x)))))
        
    # set cut if not set
    if cut is 'fixed':
        cut = max(max(df[key]), max(-df[key])) / 3

    # FIXME: cut should be either:
    # value applied to all
    # 'fixed' (default) computing a shared cutoff for all rows
    # 'adaptive' computing a cutoff that fits each row.


    # make the data frame to plot
    row_iter = df.itertuples()
    col_iterators = zip(*(horizon(row, key, cut) for row in row_iter))
    col_names = ('yp1', 'yp2', 'yp3', 'yp4', 
                 'yn1', 'yn2', 'yn3', 'yn4')
    df2 = (df.copy(deep=False)
           .assign(**dict(zip(col_names, col_iterators)))
          )

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
        chrom_subset = df.loc[df.chrom == chrom_name]
        est_chrom_len = np.max(chrom_subset.start) + width
        chrom_sizes.append(est_chrom_len)
        
    # relative width of each plot facet 
    # (using lengths of chromosomes)
    facet_widths_ratios = chrom_sizes# * nr_pop

    # make the plot
    with sns.axes_style("ticks"):

        # make the facet grid
        g = sns.FacetGrid(df2, 
                          col=chrom, 
                          row=pop,
                          # sharex=False,
                          sharex=True,
                          # margin_titles=True,
                          size=0.5, 
                          aspect=40,
                          col_order=sorted_chrom_names,
                          row_order=pop_sorting,                      
                          gridspec_kws={'hspace':0.0, 
                                        "width_ratios": facet_widths_ratios}
                         )

        # plot colors
        colours = sns.color_palette("Blues", 3) + ['black'] + \
                  sns.color_palette("Reds", 3) + ['grey']

        # first y tick
        ytic1 = round_to_1_signif(cut / 3)

        for col_name, colour in zip(col_names, colours):
            plt.setp(g.fig.texts, text="") # hack to make y facet labels align...
            # map barplots to each facet
            g.map(plt.bar, 
                  start, 
                  col_name, 
                  edgecolor = "none", 
                  width=width, 
                  color=colour)
            # no tick labels on x
            # g.set(xticklabels=[])
            #g.set_titles('{col_name}', '{row_name}')
    
        g.set_ylabels('')

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
                ax.set_yticks([ytic1, ytic1*2, ytic1*3])
#                ax.set_yticks([])
                g.set_titles('', '')
              
        # remove top and right frame
        sns.despine()

        #plt.tight_layout()

        plt.subplots_adjust(right=0.95)
        
        return g.fig

if __name__ == "__main__":
    
    class SwapContext():

        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
        
        def __enter__(self):
            self.orig = sns.plotting_context()
            sns.set_context(*self.args, **self.kwargs)
            
        def __exit__(self, type, value, traceback):
            sns.set_style(self.orig)


    n = 1500
    df1 = pd.DataFrame({'chrom': ['chr1']*2*n + ['chr1']*2*n + ['chr1']*2*n + ['chr1']*2*n,
                       'pop': ['hu']*1*n + ['ha']*1*n + ['hi']*1*n + ['ho']*1*n + ['hi']*1*n + ['baz']*1*n + ['bar']*1*n + ['foo']*1*n, 
                       'start': list(range(1*n)) * 8, 
                       'pi': list(np.sin(np.linspace(-np.pi, np.pi, 1*n))) * 8})

    df2 = pd.DataFrame({'chrom': ['chr2']*2*n + ['chr2']*2*n + ['chr2']*2*n + ['chr2']*2*n,
                       'pop': ['hu']*1*n + ['ha']*1*n + ['hi']*1*n + ['ho']*1*n + ['hi']*1*n + ['baz']*1*n + ['bar']*1*n + ['foo']*1*n, 
                       'start': list(range(1*n)) * 8, 
                       'pi': list(np.sin(np.linspace(-np.pi, np.pi, 1*n))) * 8})

    df = pd.concat([df1, df2])

    print(df.head())

    with Timer() as t:
        with SwapContext("notebook", font_scale=0.5):
            fig = horizon_plot(df, 'pi', width=1, pop='pop')#, cut=0.32)
            print('done')
            # save to file
    plt.savefig('tmp.pdf')
    #         fig.clf() # clean up memory
    plt.close(fig)  # close to allow garbage collection, also suppresses inline plot
    #         gc.collect()


    