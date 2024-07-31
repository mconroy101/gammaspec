import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from scipy.special import erf, erfc
from scipy.integrate import quad
import sys
from uncertainties import ufloat
import datetime

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

FIT_MODE = 1
SAVED = False


def polynomial_bg(x, a, b, c, *args):
    return a*x**2 + b*x + c

def step(x, x0, sigma, N, stp):
    return (N/(sigma*np.sqrt(2*np.pi)))*(stp/(1+np.exp((x-x0)/sigma)))

def gaussian(x, x0, sigma, N):
    return (N/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-x0)**2) / (2*sigma**2))    

def skewed_gaussian(x, x0, sigma, N, skw):
    return gaussian(x, x0, sigma, N)*(1+erf(skw*(x-x0)/(sigma*np.sqrt(2))))

def tail(x, x0, sigma, N, beta):
    return N * np.exp((x-x0)/beta) * erfc((x-x0)/(np.sqrt(2)*sigma) + sigma/(np.sqrt(2)*beta))

def fit_function_0(x, a, b, c, stp, sigma, x0, N, skw, *args):
    """
    Fits skewed gaussian with polynomial background + step.
    """
    
    fit = polynomial_bg(x, a, b, c) + step(x, x0, sigma, N, stp) + skewed_gaussian(x, x0, sigma, N, skw)
    for i in range(int(len(args)/3)):
        # print(args)
        x02, N2, skw2 = args[i*3:i*3+3]
        gaussian2 = (N2/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-x02)**2) / (2*sigma**2))
        sk_gaussian2 = gaussian2*(1+erf(skw2*(x-x02)/(sigma*np.sqrt(2))))
        step2 = (N2/(sigma*np.sqrt(2*np.pi)))*(stp/(1+np.exp((x-x02)/sigma)))
        fit += skewed_gaussian(x, x02, sigma, N2, skw2)
        fit += step(x, x02, sigma, N2, stp)

    return fit

def fit_function_1(x, a, b, c, stp, sigma, x0, N, beta, *args):
    """
    Fits gaussian with low E tail, polynomial background + step.
    """

    fit = polynomial_bg(x, a, b, c) + step(x, x0, sigma, N, stp) + gaussian(x, x0, sigma, N) + tail(x, x0, sigma, N, beta)

    for i in range(int(len(args)/3)):
        x02, N2, beta2 = args[i*3:i*3+3]
        fit += gaussian(x, x02, sigma, N2)
        fit += tail(x, x02, sigma, N2, beta2)
        fit += step(x, x02, sigma, N2, stp)

    return fit

def fit_function_2(x, a, b, c, sigma, x0, N, *args):
    """
    Fits gaussian with polynomial background.
    """

    fit = polynomial_bg(x, a, b, c) + gaussian(x, x0, sigma, N)
    for i in range(int(len(args)/2)):
        # print(args)
        x02, N2 = args[i*2:i*2+2]
        fit += gaussian(x, x02, sigma, N2)

    return fit

# Convert channels to energy
def calibrate(x, a, b, c):
    E = a + b*x + c*x*x
    return E

# Function that runs on mouse click on Matplotlib plot
def mouse_event(event):
    # Edit global variables
    global x
    global y

    if event.button == 3: # Right mouse button
        x.append(event.xdata) # Get x position of click into list
        y.append(event.ydata) # Get y position of click into list
        if len(x) == 1:
            # Draw a black vertical line at click position for limits
            plt.axvline(event.xdata, linestyle='--', color='k') 
            plt.title('Select upper bound')
            plt.draw()
        elif len(x) == 2:
            # Draw a black vertical line at click position for limits
            plt.axvline(event.xdata, linestyle='--', color='k') 
            plt.title('Select peak(s) position(s)')
            plt.draw()
            # Draw a red vertical line at click position for peaks
        elif len(x) > 2:
            plt.axvline(event.xdata, linestyle='--', color='r')
            plt.draw()

# Function that runs on keyboard press while matplotlib plot is open
def key_event(event):
    # Edit global variables
    global close
    global x
    global y
    # Close plot and send exit signal if escape key pressed
    if event.key == 'escape':
        plt.close()
        close = True
    # Otherwise, unless enter key is pressed, do nothing
    elif event.key != 'enter':
        pass
    # BELOW ALL RUN IF ENTER PRESSED
    elif len(x) == 0: # If no point selected
        print('Select lower limit')

    elif len(x) == 1: # If one point selected
        print('Select upper limit')

    elif len(x) < 3: # If two points selected
        print('Select at least one peak')

    elif x[1] < x[0]: # If limits are wrong way around
        print('Invalid limits, select low then high')
    else: # Close figure for fitting
        plt.close()

# Function that runs on keyboard press while another matplotlib plot is open 
def key_event_2(event):
    # Edit global variables
    global close 
    global SAVED
    if event.key == 'escape': # Close plot and send exit signal if escape key pressed
        plt.close()
        close = True    
    elif event.key == 'enter': # Close plot for further selection if enter key pressed
        plt.close()
    elif event.key == ' ':
        if SAVED == False:
            save_fit()
            SAVED = True
            print('Paramters saved to file.')
            # textlabel = plt.text(min(x_plot)-9, 1.08*max(y_plot), "SAVED PARAMETERS", size=16,
            # ha="left", va="top",
            # bbox=dict(boxstyle="round",
            #         ec=(1., 0.5, 0.5),
            #         fc=(1., 0.8, 0.8),
            #         )
            # )
            # textlabel.set_visible(True)
            # plt.title('Saved')
            # plt.draw()
        else:
            print('Fit already saved')
    else:
        pass

def calc_chi_sq(x_data, y_data, popt):
    # Calculate reduced chi square
    N_params = len(popt)
    if FIT_MODE == 0:
        residuals = (y_data[(x_data>=x[0])&(x_data<=x[1])] - fit_function_0(x_data[(x_data>=x[0])&(x_data<=x[1])], *popt))
    elif FIT_MODE == 1:
        residuals = (y_data[(x_data>=x[0])&(x_data<=x[1])] - fit_function_1(x_data[(x_data>=x[0])&(x_data<=x[1])], *popt))
    elif FIT_MODE == 2:
        residuals = (y_data[(x_data>=x[0])&(x_data<=x[1])] - fit_function_2(x_data[(x_data>=x[0])&(x_data<=x[1])], *popt))
    chisq = np.sum((residuals**2/y_data[(x_data>=x[0])&(x_data<=x[1])]))
    chisq_ndf = chisq/(len(y_data[(x_data>=x[0])&(x_data<=x[1])])-N_params)
    return chisq_ndf

def print_parameters(output_dict, chisq):
    print('------------\nFIT RESULTS:\n------------') 
    print(f'Chi-sq/ndf = {chisq:.5f}')
    print(f"Background: a = {output_dict['a']:.1uS}, b = {output_dict['b']:.1uS}, c = {output_dict['c']:.1uS}")
    if FIT_MODE == 0:
        print(f"Step = {output_dict['Step']:.1uS}\n")
        for i in range(len(x)-2):
            integral = quad(gaussian, x[0], x[1], args=(output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n + output_dict[f'Skew {i+1}'].n))[0]
            output_dict[f'Area {i+1}'] = ufloat(integral, np.sqrt(integral))
            output_dict[f'FWHM'] = 2.35*output_dict['Sigma']

            print(f"Peak {i+1}: Centroid = {output_dict[f'Centroid {i+1}']:.1uS}, Area = {output_dict[f'Area {i+1}']:.1uS}, Skew = {output_dict[f'Skew {i+1}']:.1uS}")

    elif FIT_MODE == 1:
        print(f"Step = {output_dict['Step']:.1uS}\n")
        for i in range(len(x)-2):
            integral = quad(gaussian, x[0], x[1], args=(output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n))[0] + \
                quad(tail, x[0], x[1], args=(output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n, output_dict[f'Beta {i+1}'].n))[0]
            output_dict[f'Area {i+1}'] = ufloat(integral, np.sqrt(integral))
            output_dict[f'FWHM'] = 2.35*output_dict[f'Sigma']
        
            print(f"Peak {i+1}: Centroid = {output_dict[f'Centroid {i+1}']:.1uS}, Area = {output_dict[f'Area {i+1}']:.1uS}, Beta = {output_dict[f'Beta {i+1}']:.1uS}")

    elif FIT_MODE == 2:
        for i in range(len(x)-2):
            integral = quad(gaussian, x[0], x[1], args=(output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n))[0]
            output_dict[f'Area {i+1}'] = ufloat(integral, np.sqrt(integral))
            output_dict[f'FWHM'] = 2.35*output_dict[f'Sigma']
        
            print(f"Peak {i+1}: Centroid = {output_dict[f'Centroid {i+1}']:.1uS}, Area = {output_dict[f'Area {i+1}']:.1uS}")

    output_dict['Chisq'] = chisq

    return output_dict

def save_fit():
    # Write line to file in this order:
    #   DATETIME    FILE_NAME    FIT_MODE    CENTROID    FWHM    AREA    (SKEW/BETA)
    try:
        with open('fit_results.csv', 'x') as f:
            f.write('DATETIME\tFILE_NAME\tFIT_MODE\tCENTROID\tFWHM\tAREA\tSKEW/BETA\n')
    except FileExistsError:
        pass
    with open('fit_results.csv', 'a') as f:
        # print(f.readlines())
        for i in range(len(x)-2):
            if FIT_MODE == 0:
                write_line = f"{datetime.datetime.now().strftime('%c')}\t{file_name}\t{FIT_MODE}\t{output_dict[f'Centroid {i+1}']}\t{output_dict[f'FWHM']}\t{output_dict[f'Area {i+1}']}\t{output_dict[f'Skew {i+1}']}\n"
            elif FIT_MODE == 1:
                write_line = f"{datetime.datetime.now().strftime('%c')}\t{file_name}\t{FIT_MODE}\t{output_dict[f'Centroid {i+1}']}\t{output_dict[f'FWHM']}\t{output_dict[f'Area {i+1}']}\t{output_dict[f'Beta {i+1}']}\n"
            elif FIT_MODE == 2:
                write_line = f"{datetime.datetime.now().strftime('%c')}\t{file_name}\t{FIT_MODE}\t{output_dict[f'Centroid {i+1}']}\t{output_dict[f'FWHM']}\t{output_dict[f'Area {i+1}']}\tn/a\n"
            f.write(write_line)
    

if __name__ == "__main__":
    # Display welcome message
    print('Welcome to peakfit.py by Max Conroy.\nCurrently, this program can fit single skewed gaussians with a polynomial and step function background.\n\nInstructions:\n\
        - First, use MatPlotLib to zoom in on the peak you wish to fit\n\
        - Use the right mouse button to select the lower and upper limits to fit\n\
        - Use the right mouse button to select the centre of the peak to fit \n\
        - Press enter to fit the peak. A new window of the fitted peak will be displayed.\n\
        - Press space to save the fit parameters to a file "fit_results.csv"\n\
        - Press enter again to close the fit plot and return to the peak selection\n\
        - Press escape at any time to close the program')
    
    print('\nFit Modes:\n0) Skewed Gaussian + Step + Polynomial\n1) Gaussian + Skewed Gaussian Tail + Step + Polynomial\n2) Gaussian + Step + Polynomial')
    
    try:
        file_name = str(sys.argv[1])
    except IndexError:
        print('No file name given, please enter a file name:')
        file_name = str(input("> "))

    try:
        FIT_MODE = int(sys.argv[2])
    except IndexError:
        print('No fit mode provided, please enter an integer fit mode:')
        FIT_MODE = int(input("> "))

    y_data = np.genfromtxt(file_name, skip_header=12, skip_footer=24)
    x_data = np.arange(len(y_data))

    # Boolean flag as to stop program when plot closes
    close = False
    while True:
        SAVED = False
        x = []
        y = []
        # Show full spectrum  and map mouse and keyboard events to relevant functions
        fig = plt.figure(figsize=(16,9))
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized() # Display fullscreen
        click = fig.canvas.mpl_connect('button_press_event', mouse_event)
        key = fig.canvas.mpl_connect('key_press_event', key_event)
        edges = np.array(x_data)-np.diff(x_data)[0]/2
        plt.suptitle('PEAK SELECTION')
        plt.title('Select lower bound')
        plt.stairs(y_data[:-1], edges)
        plt.show(block=True)
        # Exit program if flagged
        if close == True:
            break
        
        # DO FITTING
        if FIT_MODE == 0:

            # Initial fit parameters for a single peak [a, b, c, stp, x0, sigma, N, skew]
            p0 = [-0.08, 11, 91, 0, 2, x[2], 150000, 0.8]
            output_strings = ['a', 'b', 'c', 'Step', 'Sigma', 'Centroid 1', 'N 1', 'Skew 1']
            # Add [x0, N, skew] to p0 for each additional peak selected
            for i in range(len(x)-3):
                p0 = p0 + [x[3+i], 150000, 0.8]
                output_strings = output_strings + [f'Centroid {i+2}', f'N {i+2}', f'Skew {i+2}']
            # Perform fit of seleted peak, over selected region 
            popt, pcov = curve_fit(fit_function_0, x_data[(x_data>=x[0])&(x_data<=x[1])], y_data[(x_data>=x[0])&(x_data<=x[1])], p0=p0, sigma=np.sqrt(y_data[(x_data>=x[0])&(x_data<=x[1])]), absolute_sigma=True)

        elif FIT_MODE == 1:
            print(len(x))
            # Initial fit parameters for a single peak
            p0 = [-0.08, 11, 91, 0, 2, x[2], 150000, 0.8]
            output_strings = ['a', 'b', 'c', 'Step', 'Sigma', 'Centroid 1', 'N 1', 'Beta 1']
            # Add [x0, N, beta] to p0 for each additional peak selected
            for i in range(len(x)-3):
                p0 = p0 + [x[3+i], 150000, 0.8]
                output_strings = output_strings + [f'Centroid {i+2}', f'N {i+2}', f'Beta {i+2}']
            # Perform fit of seleted peak, over selected region 
            popt, pcov = curve_fit(fit_function_1, x_data[(x_data>=x[0])&(x_data<=x[1])], y_data[(x_data>=x[0])&(x_data<=x[1])], p0=p0, sigma=np.sqrt(y_data[(x_data>=x[0])&(x_data<=x[1])]), absolute_sigma=True)
                
        elif FIT_MODE == 2:

            # Initial fit parameters for a single peak
            p0 = [-0.08, 11, 91, 2, x[2], y[2]]
            output_strings = ['a', 'b', 'c', 'Sigma', 'Centroid 1', 'N 1']
            # Add [x0, N, skew] to p0 for each additional peak selected
            for i in range(len(x)-3):
                p0 = p0 + [x[3+i], 4000]
                output_strings = output_strings + [f'Centroid {i+2}', f'N {i+2}']
            # Perform fit of seleted peak, over selected region 
            popt, pcov = curve_fit(fit_function_2, x_data[(x_data>=x[0])&(x_data<=x[1])], y_data[(x_data>=x[0])&(x_data<=x[1])], p0=p0, sigma=np.sqrt(y_data[(x_data>=x[0])&(x_data<=x[1])]), absolute_sigma=True)
            
        # Calculate chi-sq
        chisq = calc_chi_sq(x_data, y_data, popt)
        # Calculate errors
        perr = np.sqrt(np.diag(pcov)) 

        # Print fit parameters
        output_dict = {output_strings[i]: ufloat(popt[i],perr[i]) for i in range(len(output_strings))}
        output_dict = print_parameters(output_dict, chisq)

        # Data to plot (not plotted, used to determine limits)
        x_plot = x_data[(x_data>=x[0])&(x_data<=x[1])]
        y_plot = y_data[(x_data>=x[0])&(x_data<=x[1])]
        
        # Smooth points to plot function
        x_smooth = np.linspace((x[0]), x[1], len(x_data)*10)

        fig = plt.figure(figsize=(16,9))
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        key = fig.canvas.mpl_connect('key_press_event', key_event_2)

        plt.title('Fit result')
        plt.stairs(y_data[:-1], edges, color='blue')
        plt.xlim(min(x_plot)-10, max(x_plot)+10)
        plt.ylim(-50, 1.1*max(y_plot))

        # Color map for peaks
        cmap = mpl.colormaps['summer']
        colors = cmap(np.linspace(0, 0.5, int(len(x)-2)))

        if FIT_MODE == 0:
            plt.plot(x_smooth, fit_function_0(x_smooth, *popt), color='red', label=f'Fit: $\chi^2_n = $ {chisq:.2f}')
            for i in range(len(x)-2):
                plt.plot(x_smooth, skewed_gaussian(x_smooth, output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n, output_dict[f'Skew {i+1}'].n), color=colors[i], linestyle='--', label = f'Peak {i+1}')
        elif FIT_MODE == 1:
            plt.plot(x_smooth, fit_function_1(x_smooth, *popt), color='red', label=f'Fit: $\chi^2_n = $ {chisq:.2f}')
            for i in range(len(x)-2):
                plt.plot(x_smooth, gaussian(x_smooth, output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n), color=colors[i], linestyle='--', label = f'Peak {i+1}')
                plt.plot(x_smooth, tail(x_smooth, output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n, output_dict[f'Beta {i+1}'].n), color=colors[i], linestyle=':', label = f'Tail {i+1}')
        elif FIT_MODE == 2:
            plt.plot(x_smooth, fit_function_2(x_smooth, *popt), color='red', label=f'Fit: $\chi^2_n = $ {chisq:.2f}')
            for i in range(len(x)-2):
                plt.plot(x_smooth, gaussian(x_smooth, output_dict[f'Centroid {i+1}'].n, output_dict[f'Sigma'].n, output_dict[f'N {i+1}'].n), color=colors[i], linestyle='--', label = f'Peak {i+1}')

        # Plot peak and background separately
        plt.plot(x_smooth, polynomial_bg(x_smooth, *popt), color='orange', linestyle='--', label='Background')
        plt.legend()
        plt.show()

        # Exit if flagged
        if close == True:
            break