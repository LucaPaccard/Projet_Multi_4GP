import multiprocessing
from openQCM.core.ringBuffer import RingBuffer
from openQCM.core.constants import Constants
from openQCM.common.fileStorage import FileStorage
from openQCM.common.logger import Logger as Log
from openQCM.common.switcher import Overtone_Switcher_5MHz, Overtone_Switcher_10MHz
from time import time
import serial
from serial.tools import list_ports
import numpy as np
from numpy import loadtxt
from scipy.interpolate import UnivariateSpline

from progressbar import Bar, Percentage, ProgressBar, RotatingMarker,Timer

from openQCM.util.ReadLine import ReadLine as rl

from time import sleep

from numpy import loadtxt

TAG = ""#"[Multiscan]"


class MultiscanProcess(multiprocessing.Process):

    # BASELINE CORRECTION
    def baseline_correction(self,x,y,poly_order):
        
        # Estimate Baseline with Least Squares Polynomial Fit (LSP)
        coeffs = np.polyfit(x,y,poly_order)
        # Evaluate a polynomial at specific values
        poly_fitted = np.polyval(coeffs,x) 
        return poly_fitted,coeffs    

    # BASELINE COEFFICIENTS   
    def baseline_coeffs(self):
        
        # initializations
        self.polyfitted_all = None
        self.coeffs_all = None
        self.polyfitted_all_phase = None
        self.coeffs_all_phase = None
        
        # loads Calibration (baseline correction) from file
        (self.freq_all,self.mag_all,self.phase_all) = self.load_calibration_file()
        
        # Baseline correction: input signal Amplitude (sweep all frequencies)
        (self.polyfitted_all,self.coeffs_all)=self.baseline_correction(self.freq_all,self.mag_all,8)
        self.mag_beseline_corrected_all= self.mag_all-self.polyfitted_all
        
        # Baseline correction: input signal Phase (sweep all frequencies)
        (self.polyfitted_all_phase,self.coeffs_all_phase)=self.baseline_correction(self.freq_all,self.phase_all,8)
        self.phase_beseline_corrected_all= self.phase_all-self.polyfitted_all_phase 
        return self.coeffs_all

    # SAVITZKY - GOLAY FOLTER 
    def savitzky_golay(self,y, window_size, order, deriv=0, rate=1):
        
        """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
        The Savitzky-Golay filter removes high frequency noise from data.
        It has the advantage of preserving the original shape and
        features of the signal better than other types of filtering
        approaches, such as moving averages techniques.
        Parameters
        ----------
        y : array_like, shape (N,) the values of the time history of the signal.
        window_size : int the length of the window. Must be an odd integer number.
        order : int the order of the polynomial used in the filtering.
                Must be less then `window_size` - 1.
        deriv: int the order of the derivative to compute (default = 0 means only smoothing)
        Returns
        -------
        ys : ndarray, shape (N) the smoothed signal (or it's n-th derivative).
        Notes
        -----
        The Savitzky-Golay is a type of low-pass filter, particularly
        suited for smoothing noisy data. The main idea behind this
        approach is to make for each point a least-square fit with a
        polynomial of high order over a odd-sized window centered at
        the point.
        Examples
        --------
        t = np.linspace(-4, 4, 500)
        y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
        ysg = savitzky_golay(y, window_size=31, order=4)
        import matplotlib.pyplot as plt
        plt.plot(t, y, label='Noisy signal')
        plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
        plt.plot(t, ysg, 'r', label='Filtered signal')
        plt.legend()
        plt.show()
        References
        ----------
        .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
           Data by Simplified Least Squares Procedures. Analytical
           Chemistry, 1964, 36 (8), pp 1627-1639.
        .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
           W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
           Cambridge University Press ISBN-13: 9780521880688
        """
        import numpy as np
        from math import factorial
        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order)) 
        except ValueError as msg:
            raise ValueError("WARNING: window size and order have to be of type int!")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("WARNING: window size must be a positive odd number!")
        if window_size < order + 2:
            raise TypeError("WARNING: window size is too small for the polynomials order!")
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        # pad the signal at the extremes with values taken from the signal itself
        firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve( m[::-1], y, mode='valid')


    ###########################################################################
    # Resonance Frequency, Resonance Peak, Bandwidth and Q-factor/Dissipation 
    ###########################################################################
    def parameters_finder(self,freq,signal,percent):
        
        # VER 0.1.2
        # [BUG] dissipation calculation fails when the sweep signal goes under zero level
        
        # Developed a new algorithm for the calculation of dissipation: 
		# - shift the signal up make the gain resonance curve positive 
		# - select the bandwith at 50% threshold of the maximum gain, considering only the right-side of the resonance curve
		# - multiply x2 the bandwith	
        
        # find minimum 
        f_min = np.min(signal)
        # find index of minimum 
        i_min = np.argmin(signal,axis=0)
        # rescale the signal 
        signal = signal - f_min
        
        f_max = np.max(signal)          # Find maximum
        i_max= np.argmax(signal,axis=0) # Find index of maximum
        
        # setup the index for finding the leading edge
        index_m = i_max
        # loop until the index at FWHM/others is found
        while signal[index_m] > percent*f_max:
            if index_m < 1:
               # print(TAG, 'WARNING: Left value not found')
               self._err1 = 1
               break
            index_m = index_m-1     
        #linearly interpolate between the previous values to find the value of freq at the leading edge
        m = (signal[index_m+1] - signal[index_m])/(freq[index_m+1] - freq[index_m])
        c = signal[index_m] - freq[index_m]*m
        i_leading = (percent*f_max - c)/m
        
        # setup index for finding the trailing edge
        index_M = i_max
        
        # loop until the index at FWHM/others is found
        while signal[index_M] > percent*f_max:
            if index_M >= len(signal)-1:
                print(TAG, 'WARNING: Right value not found')
                self._err2 = 1
                break
            index_M = index_M+1;
        
        # linearly interpolate between the previous values to find the value of freq at the trailing edge
        m = (signal[index_M-1] - signal[index_M])/(freq[index_M-1] - freq[index_M])
        c = signal[index_M] - freq[index_M]*m
        i_trailing = (percent*f_max - c)/m
        
        #compute the FWHM/others
        # bandwidth = abs(i_trailing - i_leading)
        
        # bandwidth = abs(i_trailing - i_leading)
        
        # VER 0.1.2
        # select the bandwith at 50% threshold of the maximum gain, considering only the right-side of the resonance curve
        # multiply x2 the bandwith	
        bandwidth = (freq[index_M] - freq[i_max])*2
        
        Qfac=freq[i_max]/bandwidth
        return i_max, f_max, bandwidth, index_m, index_M, Qfac  

    
    # ELABORATE SIGNAL 
    # -------------------------------------------------------------------------
    # TODO elaborate multi signal     
    def elaborate_multi(self, k, overtone_number, coeffs_all, readFREQ, samples, 
                  Xm, Xp, temperature, SG_window_size, Spline_points, Spline_factor, timestamp):
        
        # Number of spline points
        points = Spline_points
        # sweep counter
        self._k= k
        # current overtones number 
        self._overtone_number = overtone_number
        # evaluated polynomial coefficients
        self._coeffs_all = coeffs_all
        # frequency range, samples number
        self._readFREQ = readFREQ
        self._samples = samples
        # support vectors
        self._Xm = Xm
        self._Xp = Xp
        self._filtered_mag = np.zeros(samples)
        # save current data 
        mag   = self._Xm
        phase = self._Xp 

        # Initializations of support vectors for later storage
        self._Xm = np.linspace(0,0,self._samples)
        self._Xp = np.linspace(0,0,self._samples)
        
        # Evaluate a polynomial at specific values based on the coefficients and frequency range
        self._polyfitted = np.polyval(self._coeffs_all, self._readFREQ)
        
        # BASELINE CORRECTION ROI (raw data)
        mag_beseline_corrected = mag-self._polyfitted
        
        # FILTERING - Savitzky-Golay
        filtered_mag = self.savitzky_golay(mag_beseline_corrected, window_size = SG_window_size, order = Constants.SG_order)
        
        # peak, index e frequency of max detection baseline corrected (filtering optional)
        #self._vector_max_baseline_corrected.append(max(mag_beseline_corrected))   #Z axis (max)
        #self._index_max_baseline_corrected.append(np.argmax(mag_beseline_corrected, axis=0)) # X axis (max position)
        #h=self._index_max_baseline_corrected.append(np.argmax(mag_beseline_corrected, axis=0))
        #self._freq_max_baseline_corrected.append(readFREQ[int(h)])
        
        # FITTING/INTERPOLATING - SPLINE
        xrange = range(len(filtered_mag))
        freq_range = np.linspace(self._readFREQ[0], self._readFREQ[-1], points)
        s = UnivariateSpline(xrange, filtered_mag, s= Spline_factor)
        xs = np.linspace(0, len(filtered_mag)-1, points)
        mag_result_fit = s(xs)
        
        # PARAMETERS FINDER
# =============================================================================
#         (index_peak_fit, max_peak_fit, bandwidth_fit, index_f1_fit, index_f2_fit, Qfac_fit) = self.parameters_finder(freq_range, mag_result_fit, percent = 0.707)
# =============================================================================
        (index_peak_fit, max_peak_fit, bandwidth_fit,index_f1_fit,index_f2_fit, Qfac_fit)= self.parameters_finder(freq_range, mag_result_fit, percent=0.5)
        
        self._my_list_f[overtone_number].append( freq_range[int(index_peak_fit)] )
        self._my_list_d[overtone_number].append( (1/Qfac_fit) )
        
        #self._temperature_buffer.append(temperature)
        self._temperature_buffer_0.append(temperature)
       
        if self._k >= self._environment:
           # FREQUENCY 
           self._vec_app1 [overtone_number] = self.savitzky_golay(self._my_list_f[overtone_number].get_all(), 
                          window_size = Constants.SG_window_environment, 
                          order = Constants.SG_order_environment)
           # TODO insert a median 
           self._freq_range_mean [overtone_number] = np.average( self._vec_app1 [overtone_number] )
           
           #DISSIPATION  
           self._vec_app1d [overtone_number] = self.savitzky_golay(self._my_list_d[overtone_number].get_all(), 
                           window_size = Constants.SG_window_environment, 
                           order = Constants.SG_order_environment)
           # TODO insert a median 
           self._diss_mean [overtone_number] = np.average( self._vec_app1d [overtone_number] )
           
           # TEMPERATURE
           if overtone_number == 0:
               self._vec_app1t = self.savitzky_golay(self._temperature_buffer_0.get_all(), 
                                                     window_size = Constants.SG_window_environment, 
                                                     order = Constants.SG_order_environment)
               self._temperature_mean = np.average(self._vec_app1t)
        
# =============================================================================
#         else:
#              # TODO necessary to avoid exception when calling elaborate_multi() when k < Constants.environment
#              self._freq_range_mean [overtone_number] = 0
#              self._diss_mean [overtone_number] = 0
#              self._temperature_mean = 0
# =============================================================================
             
        # TIME EPOCH TODO 
        # ---------------------------------------------------------------------
        import datetime
        epoch = datetime.datetime(1970, 1, 1, 0, 0) #offset-naive datetime
        ts_mult = 1e6
        
        # TODO the Time is now and it is hard 
        if overtone_number == 0:
            self._my_time = (int((datetime.datetime.now() - epoch).total_seconds()*ts_mult)) #datetime.datetime.utcnow()
        # ---------------------------------------------------------------------
        
        # time array for each harmonic
        self._my_time_array[overtone_number] = (int((datetime.datetime.now() - epoch).total_seconds()*ts_mult))
        
        # TODO ADD BUFFER MEASUREMENT DATA TO THE PARSER QUEUE
        # ------------------------------------------------------
        # AMPLITUDE 
        self._parser1.add1(filtered_mag) 
        # PHASE 
        self._parser2.add2(phase)        
        
       
        # TODO just dummy 
        # Adds "fake" frequency, dissipation and temperature meaan to parser queues
        self._parser3.add3([self._my_time,0]) 
        self._parser4.add4([self._my_time,0]) 
        self._parser5.add5([self._my_time, self._temperature_mean])
        
        # add multi overtone average date to parser queue
# =============================================================================
#         self._parser_F_multi.add_F_multi( [ self._my_time, self._freq_range_mean] )
#         self._parser_D_multi.add_D_multi( [ self._my_time, self._diss_mean] )
# =============================================================================
        
        # add multi overtone frequency - dissipation and correpsonding time array to the parser queues
        self._parser_F_multi.add_F_multi( [ self._my_time_array, self._freq_range_mean] )
        self._parser_D_multi.add_D_multi( [ self._my_time_array, self._diss_mean] )
        
        # TODO single sweeep data log 
       
    def elaborate_ampli_phase_multi(self, overtone_index, poly_coeff, freq_sweep, amp_sweep, phase_sweep):
# =============================================================================
#         # init variable TODO filtering - Savitzky-Golay
#         amp_filtered = np.zeros(Constants.SAMPLES) # unused 
#         
#         freq_np = np.linspace(0,0,self._samples)
#         amp_np = np.linspace(0,0,self._samples)
#         phase_np = np.linspace(0,0,self._samples)
# =============================================================================
        
        # calculate polynomial at frequency sweep points
        amp_baseline = np.polyval(poly_coeff, freq_sweep)
        
        self._my_list_amp [overtone_index] = (amp_sweep - amp_baseline).tolist()
        self._my_list_phase [overtone_index] = phase_sweep.tolist()
        self._my_list_freq [overtone_index ] = freq_sweep.tolist()
     
        # add to new parser 
        self._parser_A_multi.add_A_multi([ self._my_list_freq, self._my_list_amp ])
        self._parser_P_multi.add_P_multi([ self._my_list_freq, self._my_list_phase ] )

    # INIT PROCESS 
    # -------------------------------------------------------------------------
    def __init__(self, parser_process):
        """
        :param parser_process: Reference to a ParserProcess instance.
        :type parser_process: ParserProcess.
        """
        multiprocessing.Process.__init__(self)
        self._exit = multiprocessing.Event()
        
        # Instantiate a ParserProcess class for each communication channels
        self._parser1 = parser_process
        self._parser2 = parser_process
        self._parser3 = parser_process
        self._parser4 = parser_process
        self._parser5 = parser_process
        # parser process for sweep info nd utility error
        self._parser6 = parser_process
        
        
        self._dummy = True
        
        self.temperature_set_old = 0
        self.cycling_time_set_old = 0
        self.P_share_set_old = 0
        self.I_share_set_old = 0
        self.D_share_set_old = 0
        
        # self.Temperature_Pid_default = [0, 0, 0, 0, 0]
        self.Temperature_Pid_default = [self.temperature_set_old, 
                                        self.cycling_time_set_old, 
                                        self.P_share_set_old, 
                                        self.I_share_set_old, 
                                        self.D_share_set_old]
        # self.temperature_set_old = loadtxt(Constants.manual_frequencies_path)
        self.Temperature_Pid_default = loadtxt(Constants.manual_frequencies_path)
        
        #DEV 
        # control temperature switch default value 
        self.ctrl_bool_pre = 0
        
        # TODO INIT PARSER PROCESS for MULTISCAN
        # Frequency 
        self._parser_F_multi = parser_process
        # Dissipation 
        self._parser_D_multi = parser_process
        # Amplitude 
        self._parser_A_multi = parser_process
        # Phase 
        self._parser_P_multi = parser_process

        # serial process 
        self._serial = serial.Serial()
        
        # self.temperature_set_old = loadtxt(Constants.manual_frequencies_path)
        
        # init ring buffer for each harmonic
        self._environment = Constants.environment
        self._frequency_buffer_0   = RingBuffer(self._environment)
        self._dissipation_buffer_0 = RingBuffer(self._environment)
        self._temperature_buffer_0 = RingBuffer(self._environment)
        self._frequency_buffer_1   = RingBuffer(self._environment)
        self._dissipation_buffer_1 = RingBuffer(self._environment)
        self._temperature_buffer_1 = RingBuffer(self._environment)
        self._frequency_buffer_2   = RingBuffer(self._environment)
        self._dissipation_buffer_2 = RingBuffer(self._environment)
        self._temperature_buffer_2 = RingBuffer(self._environment)
        self._frequency_buffer_3   = RingBuffer(self._environment)
        self._dissipation_buffer_3 = RingBuffer(self._environment)
        self._temperature_buffer_3 = RingBuffer(self._environment)
        # TODO 5M
        self._frequency_buffer_4   = RingBuffer(self._environment)
        self._dissipation_buffer_4 = RingBuffer(self._environment)
        self._temperature_buffer_4 = RingBuffer(self._environment)
        
        # init frequency dissipation list of ring buffer 
        self._my_list_f = [self._frequency_buffer_0, self._frequency_buffer_1, self._frequency_buffer_2, self._frequency_buffer_3, self._frequency_buffer_4]
        self._my_list_d = [self._dissipation_buffer_0, self._dissipation_buffer_1, self._dissipation_buffer_2, self._dissipation_buffer_3, self._dissipation_buffer_4]

        # TODO IMPORTANT chenge the init of array 
        # TODO 5M modified the number of items in the array below
        
        # init array for frequency, dissipation and temperature
        self._vec_app1 = [0,0,0,0,0]
        self._freq_range_mean = [0,0,0,0,0]
        self._vec_app1d = [0,0,0,0,0]
        self._diss_mean = [0,0,0,0,0]
        
        self._my_time = 0
        self._my_time_array = [0,0,0,0,0]
        
        # init amplitude, phase frequency sweep buffer for each harmonic
        # fundamental 
        self._amp_sweep_0 = None
        self._phase_sweep_0 = None
        self._freq_sweep_0 = None
        # 3rd overtone 
        self._amp_sweep_1 = None
        self._phase_sweep_1 = None
        self._freq_sweep_1 = None
        # 5th overtone 
        self._amp_sweep_2 = None
        self._phase_sweep_2 = None
        self._freq_sweep_2 = None
        # 7th ovetone 
        self._amp_sweep_3 = None
        self._phase_sweep_3 = None
        self._freq_sweep_3 = None
        # 9th overtone 
        self._amp_sweep_4 = None
        self._phase_sweep_4 = None
        self._freq_sweep_4 = None
        
        # init amplitude, phase and frequency list of buffer
        self._my_list_amp = [ self._amp_sweep_0 ,  self._amp_sweep_1 ,  self._amp_sweep_2,  self._amp_sweep_3,  self._amp_sweep_4 ]
        self._my_list_phase = [ self._phase_sweep_0, self._phase_sweep_1 , self._phase_sweep_2 , self._phase_sweep_3 , self._phase_sweep_4]
        self._my_list_freq = [ self._freq_sweep_0, self._freq_sweep_1 , self._freq_sweep_2 , self._freq_sweep_3 , self._freq_sweep_4]
        
        # DEBUG_0.1.1a
        # byte available at port
        self.byte_at_port = 0
        # DEBUG_0.1.1a
        # a boolean variable to check timeout at serial port and brek the for 
        self.TIME_OUT = 0
        
    # SERIAL PORT OPEN and general setting   
    # -------------------------------------------------------------------------
    def open(self, port,  speed = Constants.serial_default_overtone, 
             timeout = Constants.serial_timeout_ms, writeTimeout = Constants.serial_writetimeout_ms):
        """
        :param port: Serial port name :type port: str.
        :param speed: Overtone selected for the analysis :type speed: str.
        :param timeout: Sets current read timeout :type timeout: float (seconds).
        :param writetTimeout: Sets current write timeout :type writeTimeout: float (seconds).
        :return: True if the port is available :rtype: bool.
        """
    
        self._serial.port = port
        self._serial.baudrate = Constants.serial_default_speed 
        self._serial.stopbits = serial.STOPBITS_ONE
        self._serial.bytesize = serial.EIGHTBITS
        self._serial.timeout = timeout
        self._serial.writetimeout = writeTimeout

        #self._overtone = float(speed)
        
        # TODO check peaks mag vector here !!!
        # Loads frequencies from file
        peaks_mag = self.load_frequencies_file()
        
# =============================================================================
#         Log.i(TAG, "OPEN SERIAL PORT")
#         print("OPEN SERIAL PORT")
# =============================================================================
        
       
        # TODO SOLVED
        # ---------------------------------------------------------------------
        # GET the SELECTED FREQUENCY OVERTONE from the MAIN WINDOW  
        try:
           self._overtone = float(speed)
        except:   
           print(TAG, "Warning: wrong frequency selection, set default to {} Hz Fundamental".format(peaks_mag[0]))
           self._overtone = peaks_mag[0]
    
        # get the index of the overtones 
        self._overtone_int = None
        for i in range(len(peaks_mag)):
            if self._overtone == peaks_mag[i]:
               self._overtone_int = i
# =============================================================================
#                print("self._overtone_int ")
#                print(self._overtone_int)
# =============================================================================
               
        # Checks for correct frequency selection
        if self._overtone_int == None:
           print(TAG, "Warning: wrong frequency selection, set default to {} Hz Fundamental".format(peaks_mag[0])) 
           self._overtone_int = 0
        # ---------------------------------------------------------------------
        
        return self._is_port_available(self._serial.port)

    # RUN PROCESS
    # -------------------------------------------------------------------------
    def run(self):
        """
        The expected format is a buffer (sweep) and a new buffer as a new sweep. 
        The method parses data, converts each value to float and adds to a queue. 
        If incoming data can't be converted to float,the data will be discarded.
        """  
        
        # init error and flag parameters
        self._flag_error = 0
        self._flag_error_usb = 0
        self._err1 = 0
        self._err2 = 0
        
        # init frequency sweep and data process param
        startF = []
        stopF = []
        stepF = []
        readF = []
        sg_window_size = []
        spline_factor = []
        spline_points= []

        # get baseline coefficient 
        coeffs_all = self.baseline_coeffs()

        # init readline object 
        # line = rl(self._serial)
        
        frequencies_file_length = 0
        
        # set initial values of sweep buffer 
# =============================================================================
#         self._amp_sweep_0 = self._zerolistmaker(Constants.SAMPLES)
#         self._phase_sweep_0 = self._zerolistmaker(Constants.SAMPLES)
# =============================================================================
        for nn in Constants.overtone_dummy:
            self._my_list_amp[nn] = self._zerolistmaker(Constants.SAMPLES)
            self._my_list_phase[nn] = self._zerolistmaker(Constants.SAMPLES)
            self._my_list_freq[nn] = self._zerolistmaker(Constants.SAMPLES)
        
        # Checks if the serial port is currently connected
        if self._is_port_available(self._serial.port):
            
            # get the number of samples 
            samples = Constants.argument_default_samples 
            
            # TODO get the number of overtones in the peak frequencies file 
            frequencies_file = self.load_frequencies_file() 
            frequencies_file_length = len(frequencies_file)
            
            # Get array sweep paramaters 
            (startF, stopF, stepF, readF, sg_window_size, spline_factor, spline_points) = self.get_frequencies(samples)
        
            # Gets the state of the serial port
            if not self._serial.isOpen(): 
                
                # open the serial port
                # -------------------------------------------------------------
                self._serial.open() 
                
                # Initializes the sweep counter
                k = 0 
                print(TAG,' Buffering and processing early raw data...')
                
                # creates a timestamp
                timestamp = time()
                
                # init frequency, dissipation and temperature ring buffer 
                self._environment = Constants.environment
                self._frequency_buffer   = RingBuffer(self._environment)
                self._dissipation_buffer = RingBuffer(self._environment)
                self._temperature_buffer = RingBuffer(self._environment)
                
                # Initializes the progress bar  
                # bar = ProgressBar(widgets=[TAG,' ', Bar(marker='>'),' ',Percentage(),' ', Timer()], maxval=self._environment).start() #
                
                # ACQUISITION LOOP
                # -------------------------------------------------------------
                while not self._exit.is_set():
                    # data reset for new sweep 
                    data_mag = np.linspace(0,0,samples)   
                    data_ph  = np.linspace(0,0,samples)
                    overtone_index = 0
                    self._boolean_buffer_length = 0
                    
                    # cycle on all overtones 
                    for overtone_index in range(len(frequencies_file)):
                        
                        # DEBUG_0.1.1a
                        # print( "for loop overtone index and param= ", overtone_index, startF[overtone_index], stopF[overtone_index], int(stepF[overtone_index]) )
                       
                        try:
                            # amplitude/phase convert bit to dB/Deg parameters
                            vmax = 3.3
                            bitmax = 4096 
                            ADCtoVolt = vmax / bitmax
                            VCP = 0.9
                            
                            # DEBUG_0.1.1a
                            # get swepp time start 
                            timeStart = time()
                            
                            # -------------------------------------------------
                            # START SWEEP 
                            
                            # WRITE SWEEP COMMAND MESSAGE TO SERIAL PORT
                            # -------------------------------------------------
                            cmd = str(startF[overtone_index]) + ';' + str(stopF[overtone_index]) + ';' + str(int(stepF[overtone_index])) + '\n'
                            self._serial.write(cmd.encode())
                            
                            # DEBUG_0.1.1a
                            # added a short sleep before read serial
                            sleep(Constants.WRITE_SERIAL_WAIT)
                            
                            # Initializes buffer and strs record
                            buffer = ''
                            strs = ["" for x in range(samples + 2)]
                        
                        except:
                            print(TAG, "Info: exception serial write fail", end='\n')
                            Log.i(TAG, "Info: exception serial write fail")
                            self._flag_error_usb = 1
                            
                        if self._flag_error_usb == 0:
                            try:
                                
                                # READ SWEEP DATA AT SERIAL PORT
                                # -------------------------------------------------
                                while 1:
                                    # append string read at serial port to buffer 
                                    
                                    # DEBUG_0.1.1a
                                    self.byte_at_port = self._serial.inWaiting()
                                    buffer += self._serial.read(self.byte_at_port).decode(Constants.app_encoding)
                                    
                                    # check the time elapsed in serial read loop
                                    _time_elapsed = time() - timeStart
                                    
                                    # print ("acquiring buffer", _time_elapsed, overtone_index, "\n" )
                                    
# =============================================================================
#                                     nlines = buffer.count('\n')
# =============================================================================
                                
                                    # TODO use read line to decrease sampling time 
                                    # buffer += line.readline().decode() 
                                    
                                    # check for EOM character 
                                    if 's' in buffer:
                                        # print (_time_elapsed)
# =============================================================================
#                                         # reset input output buffer
#                                         self._serial.reset_input_buffer()
#                                         self._serial.reset_output_buffer()
#                                         sleep(0.1)
# =============================================================================
                                        # DEBUG_0.1.1a
                                        # print ("find EOM ovetone number = ", overtone_index)
                                        # DEBUG_0.1.1a
                                        #print(TAG, "EOM: byte at port = ", self.byte_at_port, end='\n')
                                        
                                        break
                                    
                                    # DEBUG_0.1.1a
                                    # insert a timeout in while acquisition loop to prevent freezing
                                    if  _time_elapsed > Constants.TIME_ELAPSED_TIMEOUT: 
                                        # DEBUG_0.1.1a
                                        print(TAG, "Info: timeout at overtone index = ", overtone_index , end='\n')
                                        self._flag_error_usb = 1
                                        self.TIME_OUT = 1
                                        # exit the while loop
                                        break   
                                    
                                if (self.TIME_OUT == 1):
                                    # DEBUG_0.1.1a
                                    # break the for ?
                                    self.TIME_OUT = 0
                                    sleep(0.5)
                                    # reset serial input/output buffer
                                    self._serial.reset_input_buffer()
                                    self._serial.reset_output_buffer()
                                    sleep(0.5)
                                    self._flag_error_usb = 1
                                    # print(TAG, "INFO: current overtone index = ", overtone_index , end='\n')
                                    
                                    # exit the for loop to prevent overtone mixing
                                    break
# =============================================================================
#                                     if nlines > Constants.argument_default_samples + 1:
#                                         print("EXCEPTION: exceed number of lines at serial port")
#                                         Log.i(TAG, "EXCEPTION: exceed number of lines at serial port")
#                                         
#                                         # reset input output buffer
#                                         self._serial.reset_input_buffer()
#                                         self._serial.reset_output_buffer()
#                                         
#                                         print("OVERTONE INDEX = ")
#                                         print (overtone_index)
#                                         
#                                         print(buffer)
#                                         
#                                         # self._flag_error_usb = 1
#                                         self._boolean_buffer_length = 1
#                                         break
# =============================================================================
                                   
     
                                # STOP SWEEP 
                                # -----------------------------------------
                                
                                if ( self._flag_error_usb == 0 ):
                                    # meaure time elapsed 
                                    # print ("Time elapsed for sweep")
                                    # print (time() - timestamp)
                                    
                                    # split each line
                                    data_raw = buffer.split('\n')
                                    length = len(data_raw)
                                    
                                    #  check the length of the serial read buffer if exceed the number of samples = 500
                                    if length > Constants.argument_default_samples + 2:
                                        print (TAG, "Info: exceed read buffer length = ", length, end='\n')
                                        self._flag_error_usb = 1
                                        data_mag = np.linspace(0,0,samples)   
                                        data_ph  = np.linspace(0,0,samples)
                                        # reset data raw
                                        data_raw = ""
                                        # reset buffer
                                        buffer = ""
                                        
                                        # reset serial input/output buffer
                                        sleep(0.5)
                                        self._serial.reset_input_buffer()
                                        self._serial.reset_output_buffer()
                                        sleep(0.5)
                                        self._flag_error_usb = 1
                                        print(TAG, "Info: current overtone index = ", overtone_index , end='\n')
                                        # break
                                        
                                    elif length < Constants.argument_default_samples + 2:
                                        # split data via semicolon ";"
                                        for i in range (length):
                                            strs[i] = data_raw[i].split(';')
                                            
                                        # TODO check the number of lines to read
                                        
                                        # converts data values to gain and phase 
                                        for i in range (length - 1):
                                            data_mag[i] = float(strs[i][0]) * ADCtoVolt / 2
                                            data_mag[i] = (data_mag[i]-VCP) / 0.03
                                            data_ph[i] = float(strs[i][1]) * ADCtoVolt / 1.5
                                            data_ph[i] = (data_ph[i]-VCP) / 0.01
                                
# =============================================================================
#                             if (overtone_index == 0):
#                                 self._amp_sweep_0 = data_mag.tolist()
#                                 self._phase_sweep_0 = data_ph.tolist()
#                                 
#                                 _sweep_path = Constants.sweep_file_path
#                                 np.savetxt( _sweep_path, np.row_stack(self._amp_sweep_0) )
# =============================================================================
                                
# =============================================================================
#                             self._my_list_amp [overtone_index] = data_mag.tolist()
#                             self._my_list_phase [overtone_index] = data_ph.tolist()
#                             self._my_list_freq [overtone_index] = readF [overtone_index]
#                             
#                             # SAVE DATA FILE TO DEBUG 
#                             if overtone_index == 0:
#                                 _sweep_path = Constants.sweep_file_path
#                                 np.savetxt( _sweep_path, np.row_stack(self._my_list_freq [overtone_index]) )
# =============================================================================

                                        # get the temperature value from the buffer 
                                        data_temp = float((strs[length - 1][0]))
# =============================================================================
#                                     
#                                         # DEBUG_0.1.1a
#                                         # ==================================================================================
#                                         # SAVE DATA FILE TO DEBUG 
#                                         _debug_path = Constants.debug_file_path
#                                         self._my_list_amp [overtone_index] = data_mag.tolist()
#                                         self._my_list_phase [overtone_index] = data_ph.tolist()
#                                         self._my_list_freq [overtone_index] = readF [overtone_index]
#                                         np.savetxt( _debug_path, np.column_stack([self._my_list_freq [overtone_index],
#                                                                                self._my_list_amp [overtone_index], 
#                                                                                self._my_list_phase [overtone_index] ]) )
#                                         # ==================================================================================
# =============================================================================
                                    
                                
                                # DEBUG_0.1.1a 
# =============================================================================
#                                     
#                                 # SET TEMPERATURE and PID PARAMETERS
#                                 # ------------------------------------------------- 
#                                 self._Temperature_PID_control()
# =============================================================================
                            
                            except:
                                 print(TAG, "Info: exception at serial port read process, overtone =", overtone_index, end='\n')
                                 Log.i(TAG, "Info: exception at serial port read process")
                                 
                                 # DEBUG_0.1.1a
                                 # reset buffer 
                                 data_mag = np.linspace(0,0,samples)   
                                 data_ph  = np.linspace(0,0,samples)
                                 
                                 # reset data raw
                                 data_raw = ""
                                 # reset buffer
                                 buffer = ""
                                 
                                 self._flag_error_usb = 1
                            
                            # DEBUG_0.1.1a
                            try: 
                                # SET TEMPERATURE and PID PARAMETERS
                                self._Temperature_PID_control()
                                
                            except:
                                print(TAG, "Info: set temperature control failed", end='\n')
                                # Log.i(TAG, "EXCEPTION: exception at serial port read process")
                                self._flag_error_usb = 1
                            
# =============================================================================
#                                  # =================================================
#                                  # WHAT HAPPEN ?
#                                  print ("BUFFER = ")
#                                  print(buffer)
#                                  print(buffer)
#                                  print("DATA RAW ")
#                                  print(data_raw)
# =============================================================================
                                 
# =============================================================================
#                         # specify handlers for different exceptions        
#                         except ValueError:
#                             print(TAG, "WARNING (ValueError): convert raw to float failed", end='\n')
#                             #Log.w(TAG, "Warning (ValueError): convert Raw to float failed")
#                         except:
#                             #if self._flag_error_usb == 1:
#                             print(TAG, "WARNING: general exception at serial port read process", end='\n')
#                             self._flag_error_usb += 1
#                             
#                             # =================================================
#                             # WHAT HAPPEN ?
#                             print ("BUFFER = ")
#                             print(buffer)
#                             print("DATA RAW ")
#                             print(data_raw)
#                             # =================================================
#                             
#                             # close and open the serial port 
#                             # close serial port
#                             self._serial.close()
#                             sleep(0.1) 
#                             # open the serial port 
#                             self._serial.open() 
# =============================================================================
                            
# =============================================================================
#                             # TODO something is occured at the serial port 
#                             # blocking the program 
#                             self._serial.flushInput()
#                             self._serial.flushOutput()
#                             self._serial.close()
#                             self._serial.open() 
# =============================================================================
                            #Log.w(TAG, "Warning (ValueError): convert Raw to float failed")
                            
                        # -----------------------------------------------------
                        # DO THE DATA_PROCESSING PROCESS HERE !!!!!!!!!!!!!!!!!
                        # TODO add new serial data to parser queue
                        #self._parser1.add1(data_mag)
                        #self._parser2.add2(data_ph)
                        # -----------------------------------------------------
                        
                        # DATA PROCESSING 
                        # -----------------------------------------------------
                        # TODO if an exception is encountered before do not elaborate 
                        
                        if self._flag_error_usb == 0:
                            try:
                                self.elaborate_multi(k, overtone_index,coeffs_all, readF[overtone_index], 
                                               samples, data_mag, data_ph, data_temp, sg_window_size[overtone_index], spline_points[overtone_index], 
                                               spline_factor[overtone_index], timestamp) 
                                
                                self.elaborate_ampli_phase_multi(overtone_index, coeffs_all, readF[overtone_index], data_mag, data_ph)
                            
                            except ValueError:
                                self._flag_error = 1
                                #if k > self._environment:
                                #   print(TAG, "WARNING (ValueError): miscalculation")
                                # self._flag_error_usb = 1
                                print(TAG, "Info: value error elaborate data", end='\n')
                                Log.i(TAG, "Info: value error elaborate data")

                                
                            except:
                                self._flag_error = 1
                                #if k > self._environment:
                                #   print(TAG, "WARNING (ValueError): miscalculation")
                                # self._flag_error_usb = 1
                                if (k > Constants.environment):
                                    print(TAG, "Info: exception general error elaborate data", end='\n')
                                    Log.i(TAG, "Info: exception general error elaborate data")
                        
# =============================================================================
#                         # DEBUG_0.1.1a    
#                         elif self._flag_error_usb == 1:
#                             # DEBUG_0.1.1a
#                             data_mag = np.linspace(0,0,samples)   
#                             data_ph  = np.linspace(0,0,samples)
#                             data_temp = 0
#                             self.elaborate_multi(k, overtone_index,coeffs_all, readF[overtone_index], 
#                                                samples, data_mag, data_ph, data_temp, sg_window_size[overtone_index], spline_points[overtone_index], 
#                                                spline_factor[overtone_index], timestamp) 
#                             self.elaborate_ampli_phase_multi(overtone_index, coeffs_all, readF[overtone_index], data_mag, data_ph)
# =============================================================================
                        
                        # add data to the "strange" parameter queue
                        self._parser6.add6([self._err1, self._err2, k, self._flag_error_usb, overtone_index])
    
                        # DEBUG_0.1.1a
# =============================================================================
#                         if (self._flag_error_usb == 1):
#                             print(" a general error occured on overtone", overtone_index)
# =============================================================================
    
    
                        # refreshes error variables at each single overtone sweep
                        self._err1 = 0
                        self._err2 = 0
                        self._flag_error_usb = 0
                        self._boolean_buffer_length = 0 
                    
                    
                    # Increases sweep counter 
                    k += 1    
# =============================================================================
#                 if k == self._environment:
#                    bar.finish()
#                    print ("bar finished")
# =============================================================================
                
                # END ACQUISITION LOOP
                # -------------------------------------------------------------
                
                self._serial.close()
    
    def get_Temperature_set_Serial(self, value_T_set):
        return value_T_set
    
    def _TempCtrl(self):
        param = loadtxt(Constants.manual_frequencies_path)
        temperature_set = param[0]
        temperature_msg = 'T' + str(int(temperature_set)) + '\n'
        self._serial.write(temperature_msg.encode())
    
    def _Temperature_PID_control(self): 
        param = loadtxt(Constants.manual_frequencies_path)
        temperature_set = param[0]
        cycling_time_set = param[1]
        P_share_set= param[2]
        I_share_set = param[3]
        D_share_set = param[4]
        _var_bool = param[5]
        _ctrl_bool = param[6]
        
        if temperature_set != self.temperature_set_old:
            temperature_msg = 'T' + str(int(temperature_set)) + '\n'
            self._serial.write(temperature_msg.encode())  
            self.temperature_set_old = temperature_set
        
        if cycling_time_set != self.cycling_time_set_old:
            cycling_time_msg = 'C' + str(int(cycling_time_set)) + '\n'
            # wait for a while before complete the communication
            sleep(0.1)
            self._serial.write(cycling_time_msg.encode())  
            sleep(0.1)
            self.cycling_time_set_old = cycling_time_set
            
        if P_share_set != self.P_share_set_old: 
            P_Share_msg = 'P' + str(int(P_share_set)) + '\n'
            # wait for a while before complete the communication
            sleep(0.1)
            self._serial.write(P_Share_msg.encode())
            sleep(0.1)
            self.P_share_set_old = P_share_set
            
        if I_share_set != self.I_share_set_old: 
            I_Share_msg = 'I' + str(int(I_share_set)) + '\n'
            # wait for a while before complete the communication
            sleep(0.1)
            self._serial.write(I_Share_msg.encode())
            sleep(0.1)
            self.I_share_set_old = I_share_set
            
        if D_share_set != self.D_share_set_old: 
            D_Share_msg = 'D' + str(int(D_share_set)) + '\n'
            # wait for a while before complete the communication
            sleep(0.1)
            self._serial.write(D_Share_msg.encode())
            sleep(0.1)
            self.D_share_set_old = D_share_set
        
        #DEV    
        # check if temperature control is enabled 
        if (_ctrl_bool != self.ctrl_bool_pre):  # value changed 
            print("DEBUG: temperature control switch I/O ", _ctrl_bool)
            # init message 
            cmd = 'X' + str(int(_ctrl_bool)) + '\n'
            # wait for a while before complete the communication
            sleep(0.1)
            # serial write message 
            self._serial.write(cmd.encode())
            sleep(0.1)
            # store new value
            self.ctrl_bool_pre = _ctrl_bool
        
        # check if temperarure value is changed 
        if  (_var_bool) == 1.0:
            temperature_msg = 'T' + str(int(temperature_set)) + '\n'
            # wait for a while before complete the communication
            sleep(0.1)
            self._serial.write(temperature_msg.encode()) 
            sleep(0.1)
            _path = Constants.manual_frequencies_path
            # reset temperature control boolean variable
            np.savetxt( _path,  np.row_stack( [param[0], param[1], param[2], param[3],param[4], 0.0, param[6]] ), fmt='%d' )
            
    
    # STOP 
    def stop(self):
        
        # close serial port
        # DEBUG 0.1.1c
        try:
            self._serial.close()
            print ("serial COM port is closed")
        except: 
            print ("WARNING: unable to close COM port ")
        
        # Signals the process to stop acquiring data.
        try:
            self._exit.set()
            print ("exit acquisition loop ")
        except:
            print ("WARNING: unable exit acquisition loop")

    # GET FREQUENCIES 
    # TODO
    def get_frequencies(self, samples):

        # TODO 

        """
        :param samples: Number of samples :type samples: int.
        :return: overtone :rtype: float.
        :return: fStep, frequency step  :rtype: float.
        :return: readFREQ, frequency range :rtype: float list.
        """
        
        startF = []
        stopF = []
        SG_window_size = []
        spline_factor = []
        fStep = []
        spline_points= []
        readFREQ = []
        
        
        # Loads frequencies from calibration file
        peaks_mag = self.load_frequencies_file()
        # get numbers of overtones stored in calibration
        peaks_mag_length = len(peaks_mag)
        
        
#        # Checks QCS type 5Mhz or 10MHz
#        # Sets start and stop frequencies for the corresponding overtone
#        
#        if (peaks_mag[0] >4e+06 and peaks_mag[0]<6e+06):
#            switch = Overtone_Switcher_5MHz(peak_frequencies = peaks_mag)
#            # 0=fundamental, 1=3th overtone and so on
#            (overtone_name,
#             overtone_value, 
#             self._startFreq,
#             self._stopFreq,
#             SG_window_size,
#             spline_factor) = switch.overtone5MHz_to_freq_range(self._overtone_int)
#            
#            print(TAG,"openQCM Device setup: @5MHz")
        
        # 10 MHz get frequency sweep param
        if (peaks_mag[0] >9e+06 and peaks_mag[0]<11e+06):
            for i in range(peaks_mag_length):
                # get multiscan sweep paramters 
                (startF_temp, 
                 stopF_temp, 
                 SG_window_size_temp, 
                 spline_factor_temp) = self.getMultiscanParameters_10Mhz(peaks_mag, i)
                
                # assign multiscan sweep param
                # TODO I do not like manage the list in this way try numpy array here 
                startF.append(startF_temp)
                stopF.append(stopF_temp)
                SG_window_size.append(SG_window_size_temp)
                spline_factor.append(spline_factor_temp)
                
                # Sets the frequency step 
                fStep_temp = (stopF_temp - startF_temp)/(samples - 1)
                # assing value 
                fStep.append(fStep_temp)
                
                # set spline points
                spline_points_temp = int( (stopF_temp - startF_temp) ) + 1
                # assing value
                spline_points.append(spline_points_temp)
                
                # set frequencies array 
                frequencies_array_temp = np.arange(samples) * (fStep_temp) + startF_temp 
                # assign value 
                readFREQ.append(frequencies_array_temp)
                
        # 5 MHz get frequency sweep param
        if (peaks_mag[0]> 4e+06 and peaks_mag[0]<6e+06):
            for i in range(peaks_mag_length):
                # get multiscan sweep paramters 
                (startF_temp, 
                 stopF_temp, 
                 SG_window_size_temp, 
                 spline_factor_temp) = self.getMultiscanParameters_5Mhz(peaks_mag, i)
                
                # assign multiscan sweep param
                # TODO I do not like manage the list in this way try numpy array here 
                startF.append(startF_temp)
                stopF.append(stopF_temp)
                SG_window_size.append(SG_window_size_temp)
                spline_factor.append(spline_factor_temp)
                
                # Sets the frequency step 
                fStep_temp = (stopF_temp - startF_temp)/(samples - 1)
                # assing value 
                fStep.append(fStep_temp)
                
                # set spline points
                spline_points_temp = int( (stopF_temp - startF_temp) ) + 1
                # assing value
                spline_points.append(spline_points_temp)
                
                # set frequencies array 
                frequencies_array_temp = np.arange(samples) * (fStep_temp) + startF_temp 
                # assign value 
                readFREQ.append(frequencies_array_temp)
        
        
        
        # Sets the frequency step 
        # fStep = (self._stopFreq-self._startFreq)/(samples-1)
        
        # Sets spline points for fitting
        # spline_points = int((self._stopFreq-self._startFreq))+1
        
        # Sets the frequency range for the corresponding overtone
        # readFREQ = np.arange(samples) * (fStep) + self._startFreq 
        
        
        return (startF, stopF, fStep, readFREQ, SG_window_size, spline_factor, spline_points)
    
    def getMultiscanParameters_10Mhz(self, peaks, overtone):
        
        if (overtone == 0):
            start = peaks[0] - Constants.L10_fundamental
            stop  = peaks[0] + Constants.R10_fundamental
            SG_win_size = Constants.SG_window_size10_fundamental
            SP_factor = Constants.Spline_factor10_fundamental
        
        elif (overtone == 1): 
            start = peaks[1] - Constants.L10_3th_overtone
            stop  = peaks[1] + Constants.R10_3th_overtone
            SG_win_size = Constants.SG_window_size10_3th_overtone
            SP_factor = Constants.Spline_factor10_3th_overtone
        
        elif (overtone == 2): 
            start = peaks[2] - Constants.L10_5th_overtone
            stop  = peaks[2] + Constants.R10_5th_overtone
            SG_win_size = Constants.SG_window_size10_5th_overtone
            SP_factor = Constants.Spline_factor10_5th_overtone
        
        return(start, stop, SG_win_size, SP_factor)
        
    def getMultiscanParameters_5Mhz(self, peaks, overtone):
        # 5 MHz fundamental 
        if (overtone == 0):
            start = peaks[0] - Constants.L5_fundamental
            stop  = peaks[0] + Constants.R5_fundamental
            SG_win_size = Constants.SG_window_size5_fundamental
            SP_factor = Constants.Spline_factor5_fundamental
        # 15 MHz 3rd overtone     
        elif (overtone == 1): 
            start = peaks[1] - Constants.L5_3th_overtone
            stop = peaks [1] + Constants.R5_3th_overtone
            SG_win_size = Constants.SG_window_size5_3th_overtone
            SP_factor = Constants.Spline_factor5_3th_overtone 
        # 25 MHz 5th overtone     
        elif (overtone == 2): 
            start = peaks[2] - Constants.L5_5th_overtone 
            stop = peaks [2] + Constants.R5_5th_overtone
            SG_win_size = Constants.SG_window_size5_5th_overtone
            SP_factor = Constants.Spline_factor5_5th_overtone
        # 35 MHz 7th overtone 
        elif (overtone == 3): 
            start = peaks[3] - Constants.L5_7th_overtone  
            stop = peaks [3] + Constants.R5_7th_overtone
            SG_win_size = Constants.SG_window_size5_7th_overtone
            SP_factor = Constants.Spline_factor5_7th_overtone
        # 45 MHz 7th overtone
        # TODO
        elif (overtone == 4): 
            start = peaks[4] - Constants.L5_9th_overtone  
            stop = peaks [4] + Constants.R5_9th_overtone
            SG_win_size = Constants.SG_window_size5_9th_overtone
            SP_factor = Constants.Spline_factor5_9th_overtone
        
        return(start, stop, SG_win_size, SP_factor) 
    
    def get_readFREQ(self, samples, overtone):
        # init temp var
        startF = []
        stopF = []
        SG_window_size = []
        spline_factor = []
        fStep = []
        spline_points= []
        
        
        # Loads frequencies from calibration file
        peaks_mag = self.load_frequencies_file()
        # get numbers of overtones stored in calibration
        # peaks_mag_length = len(peaks_mag)
        
        # 10 MHz quartz resonators
        if (peaks_mag[0] >9e+06 and peaks_mag[0]<11e+06):
            # get multiscan sweep paramters 
            (startF_temp, 
             stopF_temp, 
             SG_window_size_temp, 
             spline_factor_temp) = self.getMultiscanParameters_10Mhz(peaks_mag, overtone)

            # assign multiscan sweep param for 10 MHz quartz resonators
            # TODO I do not like manage the list in this way try numpy array here 
            startF.append(startF_temp)
            stopF.append(stopF_temp)
            SG_window_size.append(SG_window_size_temp)
            spline_factor.append(spline_factor_temp)
            
            # Sets the frequency step 
            fStep_temp = (stopF_temp - startF_temp)/(samples - 1)
            # assing value 
            fStep.append(fStep_temp)
            
            # set spline points
            spline_points_temp = int( (stopF_temp - startF_temp) ) + 1
            # assing value
            spline_points.append(spline_points_temp)
            
            # set frequencies array 
            frequencies_array_temp = np.arange(samples) * (fStep_temp) + startF_temp 
        
        # 5 MHz quartz resonators
        if (peaks_mag[0] >4e+06 and peaks_mag[0]<6e+06):
            # get multiscan sweep paramters 
            (startF_temp, 
             stopF_temp, 
             SG_window_size_temp, 
             spline_factor_temp) = self.getMultiscanParameters_5Mhz(peaks_mag, overtone)

            # assign multiscan sweep param for 5 MHz quartz resonators
            # TODO I do not like manage the list in this way try numpy array here 
            startF.append(startF_temp)
            stopF.append(stopF_temp)
            SG_window_size.append(SG_window_size_temp)
            spline_factor.append(spline_factor_temp)
            
            # Sets the frequency step 
            fStep_temp = (stopF_temp - startF_temp)/(samples - 1)
            # assing value 
            fStep.append(fStep_temp)
            
            # set spline points
            spline_points_temp = int( (stopF_temp - startF_temp) ) + 1
            # assing value
            spline_points.append(spline_points_temp)
            
            # set frequencies array 
            frequencies_array_temp = np.arange(samples) * (fStep_temp) + startF_temp 
        
        return frequencies_array_temp
        
        
    # LOAD FREQUENIES FILE 
    @staticmethod
    def load_frequencies_file():
        data  = loadtxt(Constants.cvs_peakfrequencies_path)
        peaks_mag = data[:,0]
        #peaks_phase = data[:,1] #unused at the moment
        return peaks_mag
    
    # LOAD CALIBRATION FILE 
    def load_calibration_file(self):
        # Loads Fundamental frequency and Overtones from file
        peaks_mag = self.load_frequencies_file()
        
        ############################### 
        # TODO check the damn QCM type 
        ###############################
        
        '''
        # Checks QCS type 5Mhz or 10MHz
        if len(peaks_mag) == 5:
           filename = Constants.csv_calibration_path
        elif len(peaks_mag) == 3:
           filename = Constants.csv_calibration_path10 
        '''
        
        # Checks QCS type 5Mhz or 10MHz
        if (peaks_mag[0] >4e+06 and peaks_mag[0]<6e+06):
           filename = Constants.csv_calibration_path
        elif (peaks_mag[0] >9e+06 and peaks_mag[0]<11e+06):
           filename = Constants.csv_calibration_path10 
        
        
        data  = loadtxt(filename)
        freq_all  = data[:,0]
        mag_all   = data[:,1]
        phase_all = data[:,2]
        return freq_all, mag_all, phase_all

    # GET SWEEP PARAMTERS 
    def get_sweep_parameters(): 
        # TODO develop if necessary a simple get sweep param
        print ("GET SWEEP PARAMETERS HERE")
    
    # TODO get all overtone  
    @staticmethod
    def get_speeds():
        #:return: List of the Overtones :rtype: str list.
        # Loads frequencies from  file (path: 'common\')
        data  = loadtxt(Constants.cvs_peakfrequencies_path)
        peaks_mag = data[:,0]
        reversed_peaks_mag = peaks_mag[::-1]
        return [str(v) for v in reversed_peaks_mag]   

    # SEARCH TEENSY COM PORTS     
    @staticmethod
    def get_ports(): 
        from openQCM.common.architecture import Architecture,OSType
        if Architecture.get_os() is OSType.macosx:
            import glob
            return glob.glob("/dev/tty.usbmodem*")
        elif Architecture.get_os() is OSType.linux:
            import glob
            return glob.glob("/dev/ttyACM*")
        else:
            found_ports = []
            port_connected = []
            found = False
            ports_avaiable = list(list_ports.comports())
            for port in ports_avaiable:
                if port[2].startswith("USB VID:PID=16C0:0483"):
                    found = True
                    port_connected.append(port[0])
                #else:
                #    Gets a list of the available serial ports.
                #    found_ports.append(port[0])
            if found:
               found_ports = port_connected 
            return found_ports
    
    # COM PORT AVAILABLE 
    def _is_port_available(self, port):
        """
        :param port: Port name to be verified.
        :return: True if the port is connected to the host :rtype: bool.
        """
        for p in self.get_ports():
            if p == port:
                return True
        return False
    
    def _zerolistmaker(self, n):
        listofzeros = [0] * n
        return listofzeros
