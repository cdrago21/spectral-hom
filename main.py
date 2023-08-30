import numpy as np
import matplotlib as mpl
import pylab as plt
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

        ax[1].plot((ww - self.gaussian_a.omega_c) * 1e-9, np.abs(self.gaussian_a.amplitude_freq(ww)) ** 2,
                   linewidth=2,
                   label=r'$|\phi(\omega)|^2$',
                   color='black',
                   )
        ax[1].plot((ww - self.gaussian_a.omega_c) * 1e-9, np.abs(self.gaussian_b.amplitude_freq(ww)) ** 2,
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

        plt.tight_layout()
        fig.savefig('figures/Fig1.png', dpi=300, bbox_inches='tight')
        plt.show()
