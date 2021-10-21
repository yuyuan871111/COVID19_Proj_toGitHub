import numpy as np
import matplotlib.pyplot as plt
import re

def moving_avg(data, avg_range):
    avg_list = []
    for i in range(len(data) - avg_range):
        temp_list = data[i:i+avg_range]
        avg_list.append(np.average(temp_list))
    return avg_list

def quickplotxvg(wkdir = '.', xvgfile = '.xvg', ax = None, title = 'the title', 
                 xlab = 'the x label', ylab = 'the y label', 
                 with_moving_avg = True, avg_window = int(10),
                 datacolumn = 1, ps2ns = False, sci_y = False):
    
    x, y, header = [], [], []
    with open(f"{wkdir}/{xvgfile}") as f:
        for line in f:
            cols = line.split()
            m = re.search('\d+', cols[0])
            if not m == None:
                if ps2ns == True:
                    x.append(float(cols[0])/1000)
                else:
                    x.append(float(cols[0]))
                y.append(float(cols[datacolumn]))
            else:
                header.append(cols)
    
    if ax == None:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
    else:
        ax1 = ax
    
    ax1.set_title(title)    
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    
    if sci_y :
        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    else:
        ax1.ticklabel_format(axis="y")
    ax1.plot(x,y, c='r', label='raw data')
    
    if with_moving_avg == True:
        y_avg = moving_avg(y, avg_window)
        x_shift = int(avg_window/2)
        x_avg = x[ 0+x_shift : len(y_avg)+x_shift ]
        ax1.plot(x_avg,y_avg, c='blue', label='average data')
    ax1.legend()
    
    if not ax == None:
        return ax1

def read_xvg(wkdir = '.', xvgfile = '.xvg', datacolumn = 1, 
            ps2ns = False, do_moving_avg = False, avg_window = 1000):
    x, y, header, x_avg, y_avg = [], [], [], [], []
    
    with open(f"{wkdir}/{xvgfile}") as f:
        for line in f:
            cols = line.split()
            m = re.search('\d+', cols[0])
            if not m == None:
                if ps2ns == True:
                    x.append(float(cols[0])/1000)
                else:
                    x.append(float(cols[0]))
                y.append(float(cols[datacolumn]))
            else:
                header.append(cols)
    
    # create average data if wanted
    if do_moving_avg:
        y_avg = moving_avg(y, avg_window)
        x_shift = int(avg_window/2)
        x_avg = x[ 0+x_shift : len(y_avg)+x_shift ]
    
    return x, y, header, x_avg, y_avg

def quickplot(x, y, ax = None, title = "the title", xlab = 'x label', ylab = 'y label',
            sci_y = False, with_moving_avg = False, avg_window = 1000):
    
    if ax == None:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
    else:
        ax1 = ax
    
    ax1.set_title(title)    
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    
    if sci_y :
        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    else:
        ax1.ticklabel_format(axis="y")
    ax1.plot(x,y, c='r', label='raw data')
    
    if with_moving_avg == True:
        y_avg = moving_avg(y, avg_window)
        x_shift = int(avg_window/2)
        x_avg = x[ 0+x_shift : len(y_avg)+x_shift ]
        ax1.plot(x_avg,y_avg, c='blue', label='average data')
    ax1.legend()
    
    if not ax == None:
        return ax1



def multiple_variants_plot(title='', xlab='', ylab='', ax=None, variants_name=[], 
                            x_list=[], y_list=[], x_avg_list=[], y_avg_list=[], do_moving_avg=False,
                            alpha = 1.0, colorlist = 'default', showlabel = True): 
    if ax == None:
        fig = plt.figure(figsize=[14,8])
        ax1 = fig.add_subplot(111)
    else:
        ax1 = ax

    ax1.set_title(title)    
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax1.ticklabel_format(axis="y")
    
    for index, eachname in enumerate(variants_name):
        if not showlabel:
            eachname = '_nolegend_'
        if not colorlist == 'default':
            if do_moving_avg:
                ax1.plot(x_avg_list[index],y_avg_list[index], label=eachname, color = colorlist[index], alpha = alpha)
            else:
                ax1.plot(x_list[index],y_list[index], label=eachname, color = colorlist[index], alpha = alpha)
        else:
            if do_moving_avg:
                ax1.plot(x_avg_list[index],y_avg_list[index], label=eachname, alpha = alpha)
            else:
                ax1.plot(x_list[index],y_list[index], label=eachname, alpha = alpha)
    ax1.legend()

    if not ax == None:
        return ax1