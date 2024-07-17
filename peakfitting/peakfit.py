import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.integrate import quad

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

# Sum of polynomial and step function
def background(x, a, b, c, s, x0, sigma, N, skw, *args):
    """Returns a background polynomial + step function for plotting.

    Args:
        x (arr[float]): x-coordinates to transform.
        a (float): Square parameter in quadratic polynomial.
        b (float): Linear parameter in quadratic polynomial.
        c (float): Constant parameter in quadratic polynomial.
        s (float): Step function parameter.
        x0 (float): Peak centroid parameter, for step location.
        sigma (float): Standard deviation of Gaussian peak.
        N (int): Height of Gaussian peak.
        skw (float): Skew of peak.

    Returns:
        arr[float]: y-coordinates y = f(x)
    """
    # Define quadratic polynomial background
    polynomial = a*x**2 + b*x + c   
    # Define smoothed step function for first peak
    step = (N/(sigma*np.sqrt(2*np.pi)))*(s/(1+np.exp((x-x0)/sigma)))
    # Add step to polynomial
    fit = polynomial + step
    # Check if other peaks are present
    for i in range(int(len(args)/3)):
        # If so, get required arguments passed to function
        x02, N2, skw2 = args[i*3:i*3    +3]
        # Define the step for the extra peak
        step2 = (N2/(sigma*np.sqrt(2*np.pi)))*(s/(1+np.exp((x-x02)/sigma)))
        # Add new step to total background
        fit += step2

    return fit

# Skewed gaussian
def multiple_peaks(x, a, b, c, stp, x0, sigma, N, skw, *args):
    peaks = []
    info = []
    info.append([x0, sigma, N, skw])
    gaussian = (N/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-x0)**2) / (2*sigma**2))
    sk_gaussian = gaussian*(1+erf(skw*(x-x0)/(sigma*np.sqrt(2))))
    peaks.append(sk_gaussian)
    for i in range(int(len(args)/3)):
        # print(args)
        x02, N2, skw2 = args[i*3:i*3+3]
        gaussian2 = (N2/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-x02)**2) / (2*sigma**2))
        sk_gaussian2 = gaussian2*(1+erf(skw2*(x-x02)/(sigma*np.sqrt(2))))
        info.append([x02, sigma, N2, skw2])
        peaks.append(sk_gaussian2)

        # peaks.append(sk_gaussian+sk_gaussian2)
    
    return peaks

def single_peak(x, x0, sigma, N, skw):

    gaussian = (N/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-x0)**2) / (2*sigma**2))
    sk_gaussian = gaussian*(1+erf(skw*(x-x0)/(sigma*np.sqrt(2))))

    return sk_gaussian


def fit_function(x, a, b, c, stp, x0, sigma, N, skw, *args):
    """_summary_

    Args:
        x (_type_): _description_
        a (_type_): _description_
        b (_type_): _description_
        c (_type_): _description_
        stp (_type_): _description_
        x0 (_type_): _description_
        sigma (_type_): _description_
        N (_type_): _description_
        skw (_type_): _description_
        *args ([[x0, N, skw], ...]): Parameters for next n peaks. Sigma and step stay the same.

    Returns:
        arr[float]: _description_
    """
    polynomial = a*x**2 + b*x + c
    gaussian = (N/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-x0)**2) / (2*sigma**2))
    sk_gaussian = gaussian*(1+erf(skw*(x-x0)/(sigma*np.sqrt(2))))
    # tail = (N/(sigma*np.sqrt(2*np.pi)))*np.exp(-((skw**2)/(2*sigma**2)) + skw*((x-x0)/sigma**2))
    step = (N/(sigma*np.sqrt(2*np.pi)))*(stp/(1+np.exp((x-x0)/sigma)))
    fit = polynomial + step + sk_gaussian
    for i in range(int(len(args)/3)):
        # print(args)
        x02, N2, skw2 = args[i*3:i*3+3]
        gaussian2 = (N2/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-x02)**2) / (2*sigma**2))
        sk_gaussian2 = gaussian2*(1+erf(skw2*(x-x02)/(sigma*np.sqrt(2))))
        step2 = (N2/(sigma*np.sqrt(2*np.pi)))*(stp/(1+np.exp((x-x02)/sigma)))
        fit += sk_gaussian2
        fit += step2

    return fit

# Convert channels to energy
def calibrate(x, a, b, c):
    E = a + b*x + c*x*x
    return E

# Function that runs on mouse click on Matplotlib plot
def mouse_event(event):

    global x
    global y

    if event.button == 3: # Right mouse button
        # print('x: {} and y: {}'.format(event.xdata, event.ydata))
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
    global close 
    if event.key == 'escape': # Close plot and send exit signal if escape key pressed
        plt.close()
        close = True    
    elif event.key == 'enter': # Close plot for further selection if enter key pressed
        plt.close() 

# Read data first
# REPLACE THIS WITH USER INPUT
y_data = np.genfromtxt('eu152_mid.Spe', skip_header=12, skip_footer=24)
x_data = np.arange(len(y_data))

# Boolean flag as to stop program when plot closes
close = False

# Display welcome message
print('Welcome to peakfit.py by Max Conroy.\nCurrently, this program can fit single skewed gaussians with a polynomial and step function background.\n\nInstructions:\n\
       - First, use MatPlotLib to zoom in on the peak you wish to fit\n\
       - Use the right mouse button to select the lower and upper limits to fit\n\
       - Use the right mouse button to select the centre of the peak to fit (can select multiple, only first will fit)\n\
       - Press enter to fit the peak. A new window of the fitted peak will be displayed.\n\
       - Press enter again to return to the peak selection\n\
       - Press escape at any time to close the program')

# Allow user to continue to select peaks and fit them
while True:
    x = []
    y = []

    # Show full spectrum and map mouse and keyboard events to relevant functions
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
    
    # Initial fit parameters for a single peak
    p0 = [-0.08, 11, 91, 0, x[2], 2, 150000, 0.8]

    # Add [x0, N, skew] to p0 for each additional peak selected
    for i in range(len(x)-3):
        # print(i)
        p0 = p0 + [x[3+i], 150000, 0.8]


    # try:
    # Perform fit of seleted peak, over selected region 
    popt, pcov = curve_fit(fit_function, x_data[(x_data>=x[0])&(x_data<=x[1])], y_data[(x_data>=x[0])&(x_data<=x[1])], p0=p0, sigma=np.sqrt(y_data[(x_data>=x[0])&(x_data<=x[1])]), absolute_sigma=True)
    # except:
    #     print('Fit failed')
    # Calculate area under peak by integrating # (a,b,c,s,x0,sigma,N,gamma)
    # bg_area, bg_area_err = quad(background, x[0], x[1], args=(popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7]))
    # peak_area, peak_area_err = (quad(peak, x[0], x[1], args=(popt[4], popt[5], popt[6], popt[7])))
    
    # Output relevant parameters 
    print('----------\nFIT RESULTS:\n----------')   
    # Calculate and output reduced chi square
    N_params = len(popt)
    residuals = (y_data[(x_data>=x[0])&(x_data<=x[1])] - fit_function(x_data[(x_data>=x[0])&(x_data<=x[1])], *popt))
    chisq = np.sum((residuals**2/y_data[(x_data>=x[0])&(x_data<=x[1])]))
    chisq_ndf = chisq/(len(y_data[(x_data>=x[0])&(x_data<=x[1])])-N_params)
    print(f'Chi-sq/ndf = {chisq_ndf:.5f}')
    print(f'Background: a = {popt[0]:5f}, b = {popt[1]:5f}, c = {popt[2]:5f}, Step: {popt[3]:5f}\n')
    peak_info = []
    peak_info_errors = []
    reduced_popt = popt[3:]
    perr = np.sqrt(np.diag(pcov))
    reduced_perr = perr[3:]
    print(f'Peak 1:\nCentroid: {reduced_popt[1]:.2f} +- {reduced_perr[1]:.2f}')
    print(f'Sigma: {reduced_popt[2]:.2f} +- {reduced_perr[2]:.2f}')
    peak_area, peak_area_err = (quad(single_peak, x[0], x[1], args=(reduced_popt[1], reduced_popt[2], reduced_popt[3], reduced_popt[4])))
    print(f'Area: {int(peak_area)} +- {int(np.sqrt(peak_area) )}\n')
    # [a, b, c, stp, x0, sigma, n, skw, x0, n, skw, ...]
    extra_popt = popt[8:]
    extra_perr = perr[8:]
    print(extra_popt)
    for i in range(int(len(x)-3)):
        print(f'Peak {i+2}:\nCentroid: {extra_popt[3*i]:.2f} +- {extra_perr[3*i]:.2f}')
        print(f'Sigma: {reduced_popt[2]:.2f} +- {reduced_perr[2]:.2f}')
        peak_area, peak_area_err = (quad(single_peak, x[0], x[1], args=(extra_popt[3*i], reduced_popt[2], extra_popt[3*i+1], extra_popt[3*i+2])))
        print(f'Area: {int(peak_area)} +- {int(np.sqrt(peak_area) )}\n')


        # peak_info.append(reduced_popt[i*3:i*3+3])
        # peak_info_errors.append(reduced_perr[i*3:i*3+3])
        # print(f'For peak {i+1}:')
        # print(f'Centroid: {peak_info[i][1]:2f} +- {peak_info_errors[i][1]:2f}\nEnergy = {calibrate(peak_info[i][1], 1.08285e+00, 4.15692e-01, 1.10344e-07):2f} keV')
        # print(f'Sigma: {peak_info[i][2]:2f} +- {peak_info_errors[i][2]:2f}\nSkew: {peak_info[i][4]:5f}')
        # peak_area, peak_area_err = (quad(single_peak, x[0], x[1], args=(peak_info[i][1], peak_info[i][2], peak_info[i][3], peak_info[i][4])))
        # print(f'Area: {int(peak_area)} +- {int(np.sqrt(peak_area) )}\n')

    

    # Data to plot (not plotted, used to determine limits)
    x_plot = x_data[(x_data>=x[0])&(x_data<=x[1])]
    y_plot = y_data[(x_data>=x[0])&(x_data<=x[1])]

    np.savetxt('eu_peak.txt', (np.vstack((x_plot, y_plot)).T))
    
    # Smooth points to plot function
    x_smooth = np.linspace((x[0]), x[1], len(x_data)*10)

    fig = plt.figure(figsize=(16,9))
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.title('Fit result')
    plt.stairs(y_data[:-1], edges, color='blue')
    plt.xlim(min(x_plot)-10, max(x_plot)+10)
    plt.ylim(-50, max(y_plot)+0.1*max(y_plot))
    # Plot full function
    plt.plot(x_smooth, fit_function(x_smooth, *popt), color='red', label=f'Fit - $\chi^2_n = $ {chisq_ndf:.2f}')
    # Plot peak and background separately
    plt.plot(x_smooth, background(x_smooth, *popt), color='orange', linestyle='--', label='Polynomial + step background')
    # plt.plot(x_smooth, peak(x_smooth, *popt)[-1], color='k', linestyle=':', label = 'Skewed Gaussian')
    cmap = mpl.colormaps['summer']
    colors = cmap(np.linspace(0, 0.5, int(len(x)-2)))
    for i in range(int(len(x)-2)):
        plt.plot(x_smooth, multiple_peaks(x_smooth, *popt)[i], color=colors[i], linestyle='--', label = f'Peak {i+1}')
    # Connect keyboard to relevant functon
    key = fig.canvas.mpl_connect('key_press_event', key_event_2)
    plt.legend()
    plt.show()
    # Exit if flagged
    if close == True:
        break

    # print('\nWould you like to save these fit parameters? (y or n)')
    # save = input('> ')
    # if save == 'y':
    #     for i in range(int(len(reduced_popt)/5)):
    #         save_data = [peak_info[i][1], peak_info_errors[i][1]]
    #     np.savetxt('save_data.txt', save_data, delimiter=', ')


