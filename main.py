import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.constants import speed_of_light

twopi = 2 * np.pi

mpl.rcParams['font.family'] = 'serif'
cmfont = mpl.font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
mpl.rcParams['font.serif'] = cmfont.get_name()
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['axes.unicode_minus'] = False
fs = 16  # font size
ts = 12  # tick size


class GaussianPulse:
    """class for a Gaussian pulse"""

    def __init__(self, param_dict):
        """Constructor for the Gaussian pulse class.
        Args:
            param_dict (dict): dictionary containing the necessary parameters to perform a calculation:
                - lambda_c (float): central wavelength of the gaussian pulse (nm)
                - FWHM (float): Full Width at Half Maximum of the intensity of the gaussian pulse in time (ps)
        """
        self.lambda_c = param_dict["lambda_c"]
        self.FWHM = param_dict["FWHM"]

        self.omega_c = None
        self.sigma = None
        self.freq_range = None
        self.time_range = None

        self.calculate_derived_parameters()

    def calculate_derived_parameters(self):
        # set center frequency (rad/s)
        self.omega_c = twopi * speed_of_light * 1e9 / self.lambda_c

        # set standard deviation (Hz)
        self.sigma = 2 * np.sqrt(np.log(2)) / (self.FWHM * 1e-12)

        self.freq_range = np.arange(self.omega_c - 3.5 * self.sigma, self.omega_c + 3.5 * self.sigma,
                                    7 * self.sigma / 100)
        self.time_range = np.arange(- 3.5 / self.sigma, 3.5 / self.sigma, 7 / self.sigma / 100)

    def amplitude_time(self, time):
        """Amplitude of the gaussian pulse.

        Args:
            time (float): time (s)
        Returns:
            array(len(t)): amplitude of the gaussian pulse in time
        """
        return (self.sigma ** 2 / np.pi) ** 0.25 * np.exp(-self.sigma ** 2 * time ** 2 / 2 - 1j * self.omega_c * time)

    def amplitude_freq(self, omega):
        """Amplitude of the gaussian pulse.

        Args:
            omega (float): frequency (rad/s)
        Returns:
            array(len(omega)): amplitude of the gaussian pulse in frequency
        """
        return (1 / (np.pi * self.sigma ** 2)) ** 0.25 * np.exp(-(omega - self.omega_c) ** 2 / (2 * self.sigma ** 2))


class IndependentGaussianCoincidence:
    """class for the coincidence probability from independent photons with Gaussian amplitudes"""

    def __init__(self, param_dict):
        """Constructor for the Gaussian pulse class.

        Args:
            param_dict (dict): dictionary containing the necessary parameters to perform a
                calculation:
                - gaussian_a (class): first Gaussian pulse class
                - gaussian_b (Class): second Gaussian pulse class
                - delay (float): delay between Gaussian's
        """
        self.gaussian_a = param_dict["gaussian_a"]
        self.gaussian_b = param_dict["gaussian_b"]

        self.FWHM_max = np.maximum(self.gaussian_a.FWHM, self.gaussian_b.FWHM)

        self.delay_range = np.arange(-2.5 * self.FWHM_max, 2.5 * self.FWHM_max, 5 * self.FWHM_max / 100) * 1e-12

    def coincidence_probability(self, delay):
        temporal_exp = np.exp(-self.gaussian_a.sigma ** 2 * self.gaussian_b.sigma ** 2 * delay ** 2 / (
                self.gaussian_a.sigma ** 2 + self.gaussian_b.sigma ** 2))
        spectral_exp = np.exp(-(self.gaussian_a.omega_c - self.gaussian_b.omega_c) ** 2 / (
                self.gaussian_a.sigma ** 2 + self.gaussian_b.sigma ** 2))
        pre_factor = self.gaussian_a.sigma * self.gaussian_b.sigma / (
                self.gaussian_a.sigma ** 2 + self.gaussian_b.sigma ** 2)
        return 1 / 2 - pre_factor * temporal_exp * spectral_exp

    def plot_coincidence(self):
        if self.gaussian_a.FWHM > self.gaussian_b.FWHM:
            tt = self.gaussian_a.time_range
            ww = self.gaussian_b.freq_range
        else:
            tt = self.gaussian_b.time_range
            ww = self.gaussian_a.freq_range

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
        ax[0].plot(tt * 1e12, np.abs(self.gaussian_a.amplitude_time(tt)) ** 2,
                   linewidth=2,
                   label=r'$|\phi(t)|^2$',
                   color='black',
                   )
        ax[0].plot(tt * 1e12, np.abs(self.gaussian_b.amplitude_time(tt)) ** 2,
                   linewidth=2,
                   label=r'$|\varphi(t)|^2$',
                   color='red',
                   linestyle='--',
                   )
        ax[0].set_xlabel("t (ps)", fontsize=fs)
        ax[0].set_ylabel("Temporal intensity", fontsize=fs)
        ax[0].tick_params(axis='both', labelsize=ts)
        ax[0].legend()

        ax[1].plot((ww - self.gaussian_a.omega_c) * 1e-9 / twopi, np.abs(self.gaussian_a.amplitude_freq(ww)) ** 2,
                   linewidth=2,
                   label=r'$|\phi(\omega)|^2$',
                   color='black',
                   )
        ax[1].plot((ww - self.gaussian_a.omega_c) * 1e-9 / twopi, np.abs(self.gaussian_b.amplitude_freq(ww)) ** 2,
                   linewidth=2,
                   label=r'$|\varphi(\omega)|^2$',
                   color='red',
                   linestyle='--',
                   )
        ax[1].set_xlabel(r"$(\omega - \overline{\omega}_a)/2\pi$ (GHz)", fontsize=fs)
        ax[1].set_ylabel("Spectral intensity", fontsize=fs)
        ax[1].tick_params(axis='both', labelsize=ts)
        ax[1].legend()

        ax[2].plot(self.delay_range * 1e12, self.coincidence_probability(self.delay_range),
                   linewidth=2,
                   color='black')
        ax[2].set_xlabel(r"$\tau$ (ps)", fontsize=fs)
        ax[2].set_ylabel(r"$p_\mathrm{II}^\mathrm{Gauss}$", fontsize=fs)
        ax[2].tick_params(axis='both', labelsize=ts)
        ax[2].axhline(0, color='black', alpha=0.5)
        ax[2].axvline(0, color='black', alpha=0.5)

        plt.tight_layout()
        fig.savefig('figures/fig2.png', dpi=300, bbox_inches='tight')
        plt.show()






class DoubleGaussianJSA:
    """class for a double Gaussian joint-spectral amplitude"""

    def __init__(self, param_dict):
        """Constructor for the double Gaussian joint-spectral amplitude class.
        Args:
            param_dict (dict): dictionary containing the necessary parameters to perform a calculation:
                - lambda_c (float): central wavelength of the photon pairs (nm)
                - pulse_duration (float): effective pulse duration set by T_p  (ps)
                - coherent_time (float): coherence  time set by T_c  (ps)
        """
        self.lambda_c = param_dict["lambda_c"]
        self.pulse_duration = param_dict["pulse_duration"]
        self.coherence_time = param_dict["coherence_time"]

        self.omega_c = None
        self.T_p, self.T_c = None, None
        self.freq_range = None
        self.time_range = None

        self.calculate_derived_parameters()

    def calculate_derived_parameters(self):
        # set center frequency (rad/s)
        self.omega_c = twopi * speed_of_light * 1e9 / self.lambda_c

        # set T_p and T_c (s)
        self.T_p = self.pulse_duration * 1e-12
        self.T_c = self.coherence_time * 1e-12

        # set frequency range (rad/s)
        self.freq_range = np.arange(self.omega_c - 2 / self.T_c, self.omega_c + 2 / self.T_c,
                                    4 / self.T_c / 100)

    def joint_spectral_amplitude(self, omega_1, omega_2):
        """Amplitude of the gaussian pulse.

        Args:
            omega_1 (float): frequency (Hz)
            omega_2 (float): frequency (Hz)
        Returns:
            array((len(omega_1), len(omega_2)): joint spectral amplitude of the double Gaussian
        """
        pre_factor = np.sqrt(self.T_p * self.T_c / np.pi)
        left_term = np.exp(-(omega_1 - omega_2) ** 2 * self.T_c ** 2 / 4)
        right_term = np.exp(-(omega_1 + omega_2 - 2 * self.omega_c) ** 2 * self.T_p ** 2 / 4)

        return pre_factor * left_term * right_term


class DoubleGaussianCoincidence:
    """class for the coincidence probability from independent photons with Gaussian amplitudes"""

    def __init__(self, param_dict):
        """Constructor for the double Gaussian coincidence class.

        Args:
            param_dict (dict): dictionary containing the necessary parameters to perform a
                calculation:
                - double_gaussian (class): double Gaussian class
                - delay_range (float): delay range
        """
        self.double_gaussian = param_dict["double_gaussian"]

        # set delay range (s)
        self.delay_range = np.arange(-2.5 * self.double_gaussian.T_c, 2.5 * self.double_gaussian.T_c,
                                     5 * self.double_gaussian.T_c / 100)

    def coincidence_probability(self, delay):
        return 1 / 2 - 1 / 2 * np.exp(- delay ** 2 / 2 / self.double_gaussian.T_c ** 2)

    def plot_coincidence(self):
        ww = self.double_gaussian.freq_range
        W1, W2 = np.meshgrid(ww, ww)
        JSA = self.double_gaussian.joint_spectral_amplitude(W1, W2)

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(4, 8))
        im0 = ax[0].pcolormesh((ww - self.double_gaussian.freq_range) * 1e-9 / twopi,
                               (ww - self.double_gaussian.freq_range) * 1e-9 / twopi,
                                np.abs(JSA) ** 2, shading='auto') #, norm=mpl.colors.LogNorm())
        ax[0].set_xlabel(r"$(\omega_1 - \overline{\omega})/2\pi$ (GHz)", fontsize=fs)
        ax[0].set_ylabel(r"$(\omega_2 - \overline{\omega})/2\pi$ (GHz)", fontsize=fs)
        ax[0].set_title(r"Joint spectral intensity", fontsize=fs)
        ax[0].tick_params(axis='both', labelsize=ts)

        plt.colorbar(im0, cax=make_axes_locatable(ax[0]).append_axes("right", size="5%", pad=0.05))

        ax[1].plot(self.delay_range * 1e-12, self.coincidence_probability(self.delay_range),
                   linewidth=2,
                   color='black')
        ax[1].set_xlabel(r"$\tau$ (ps)", fontsize=fs)
        ax[1].set_ylabel(r"$p_\mathrm{II}^\mathrm{Gauss}$", fontsize=fs)
        ax[1].tick_params(axis='both', labelsize=ts)
        ax[1].axhline(0, color='black', alpha=0.5)
        ax[1].axvline(0, color='black', alpha=0.5)

        plt.tight_layout()
        fig.savefig('figures/fig3.png', dpi=300, bbox_inches='tight')
        plt.show()
