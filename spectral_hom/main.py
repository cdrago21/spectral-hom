import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.constants import speed_of_light

twopi = 2 * np.pi

mpl.rcParams['font.family'] = 'serif'
cmfont = mpl.font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
mpl.rcParams['font.serif'] = cmfont.get_name()
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['axes.unicode_minus'] = False
fs = 14  # font size
ts = 11  # tick size

cmap = 'GnBu'
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

    def plot_coincidence(self, figname):
        if self.gaussian_a.FWHM > self.gaussian_b.FWHM:
            tt = self.gaussian_a.time_range
            ww = self.gaussian_b.freq_range
        else:
            tt = self.gaussian_b.time_range
            ww = self.gaussian_a.freq_range

        omega_c = self.gaussian_a.omega_c
        t_lim = min(-np.min(tt), np.max(tt)) * 1e12
        t_ticks = np.round(np.linspace(-t_lim, t_lim, 7), 0)

        f_lim = max(-np.min(ww - omega_c), np.max(ww - omega_c)) * 1e-9 / twopi
        f_ticks = np.round(np.linspace(-f_lim, f_lim, 7), 0)

        tau_lim = min(-np.min(self.delay_range), np.max(self.delay_range)) * 1e12
        tau_ticks = np.round(np.linspace(-tau_lim, tau_lim, 7), 0)

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 3))
        ax[0].plot(tt * 1e12, np.abs(self.gaussian_a.amplitude_time(tt)) / np.max(np.abs(self.gaussian_a.amplitude_time(tt))),
                   linewidth=2,
                   label=r'$|\bar\phi_a(t)|$',
                   color='black',
                   )
        ax[0].plot(tt * 1e12, np.abs(self.gaussian_b.amplitude_time(tt)) / np.max(np.abs(self.gaussian_b.amplitude_time(tt))),
                   linewidth=2,
                   label=r'$|\bar\phi_b(t)|$',
                   color='red',
                   linestyle='--',
                   )
        ax[0].set_xlabel("$t$ (ps)", fontsize=fs)
        ax[0].set_ylabel(r"$|$Temporal amplitude$|$", fontsize=fs)
        ax[0].tick_params(axis='both', labelsize=ts)
        ax[0].text(0.05, 0.9, f'$\mathrm{{FWHM}}_a$={self.gaussian_a.FWHM} ps', fontsize=12, transform=ax[0].transAxes,
                   bbox = dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1'))
        ax[0].text(0.05, 0.8, f'$\mathrm{{FWHM}}_b$={self.gaussian_b.FWHM} ps', fontsize=12, transform=ax[0].transAxes,
                   bbox = dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1'))
        ax[0].set_xticks(t_ticks)
        ax[0].grid('on', alpha = 0.5)
        ax[0].legend(loc = 'upper right', prop={'size': 12})

        ax[1].plot((ww - self.gaussian_a.omega_c) * 1e-9 / twopi,
                   np.abs(self.gaussian_a.amplitude_freq(ww)) / np.max(np.abs(self.gaussian_a.amplitude_freq(ww))),
                   linewidth=2,
                   label=r'$|\phi_a(\omega)|$',
                   color='black',
                   )
        ax[1].plot((ww - self.gaussian_a.omega_c) * 1e-9 / twopi,
                   np.abs(self.gaussian_b.amplitude_freq(ww)) / np.max(np.abs(self.gaussian_b.amplitude_freq(ww))),
                   linewidth=2,
                   label=r'$|\phi_b(\omega)|$',
                   color='red',
                   linestyle='--',
                   )
        ax[1].set_xlabel(r"$(\omega - \overline{\omega}_a)/2\pi$ (GHz)", fontsize=fs)
        ax[1].set_ylabel(r"$|$Spectral amplitude$|$", fontsize=fs)
        ax[1].tick_params(axis='both', labelsize=ts)
        D = (self.gaussian_a.omega_c - self.gaussian_b.omega_c) * 1e-9 / twopi
        ax[1].text(0.05, 0.9, f"$\Delta={D:.2f}$ GHz", fontsize=12, transform=ax[1].transAxes,
                   bbox = dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1'))
        ax[1].set_xticks(f_ticks)
        ax[1].set_xlim(min(ww - self.gaussian_a.omega_c) * 1e-9 / twopi,
                   max(ww - self.gaussian_a.omega_c) * 1e-9 / twopi)
        ax[1].grid('on', alpha = 0.5)
        ax[1].legend(loc = 'upper right', prop={'size': 12})

        ax[2].plot(self.delay_range * 1e12, self.coincidence_probability(self.delay_range),
               linewidth=2,
               color='black')
        ax[2].set_xlabel(r"$\tau$ (ps)", fontsize=fs)
        ax[2].set_ylabel(r"$p_\mathrm{G}$", fontsize=fs)
        ax[2].tick_params(axis='both', labelsize=ts)
        ax[2].axhline(0, color='black', alpha=0.5)
        ax[2].axvline(0, color='black', alpha=0.5)
        ax[2].set_xticks(tau_ticks)
        ax[2].grid('on', alpha = 0.5)

        plt.tight_layout()
        fig.savefig(f"figures/{figname}.pdf", dpi=300, bbox_inches='tight')
        plt.show()


class JointGaussianJSA:
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
                                    4 / self.T_c / 400)

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


class JointGaussianCoincidence:
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
        self.delay_range = np.arange(-5 * self.double_gaussian.T_c, 5 * self.double_gaussian.T_c,
                                     10 * self.double_gaussian.T_c / 100)

    def coincidence_probability(self, delay):
        return 1 / 2 - 1 / 2 * np.exp(- delay ** 2 / 2 / self.double_gaussian.T_c ** 2)

    def plot_coincidence(self, figname1, figname2):
        ww = self.double_gaussian.freq_range
        W1, W2 = np.meshgrid(ww, ww)
        JSA = self.double_gaussian.joint_spectral_amplitude(W1, W2)

        ticks = np.linspace(int(np.min((ww - self.double_gaussian.omega_c) * 1e-9 / twopi)),
                            int(np.max((ww - self.double_gaussian.omega_c) * 1e-9 / twopi)),
                            5)

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5, 4))
        im0 = ax.pcolormesh((W1 - self.double_gaussian.omega_c) * 1e-9 / twopi,
                               (W2 - self.double_gaussian.omega_c) * 1e-9 / twopi,
                               np.abs(JSA) / np.max(np.abs(JSA)),
                               shading='gouraud',
                               cmap=cmap)
        ax.set_xlabel(r"$(\omega_1 - \overline{\omega})/2\pi$ (GHz)", fontsize=fs)
        ax.set_ylabel(r"$(\omega_2 - \overline{\omega})/2\pi$ (GHz)", fontsize=fs)
        ax.set_title(r"$|$Joint spectral amplitude$|$", fontsize=fs)
        ax.tick_params(axis='both', labelsize=ts)
        ax.set_yticks(ticks)
        ax.set_xticks(ticks)
        plt.colorbar(im0, cax=make_axes_locatable(ax).append_axes("right", size="5%", pad=0.025))
        plt.tight_layout()
        fig.savefig(f'figures/{figname1}.pdf', dpi=300, bbox_inches='tight')
        plt.show()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5, 4))
        ax.plot(self.delay_range * 1e12, self.coincidence_probability(self.delay_range),
                   linewidth=2,
                   color='black')
        ax.set_xlabel(r"$\tau$ (ps)", fontsize=fs)
        ax.set_ylabel(r"$p_\mathrm{JG}$", fontsize=fs)
        ax.tick_params(axis='both', labelsize=ts)
        ax.axhline(0, color='black', alpha=0.5)
        ax.axvline(0, color='black', alpha=0.5)
        ax.grid('on', alpha = 0.5)
        plt.tight_layout()
        fig.savefig(f'figures/{figname2}.pdf', dpi=300, bbox_inches='tight')
        plt.show()


def general_coincidence_pre(frequency1_range, dfreq1, frequency2_range, dfreq2, joint_spectral_amplitude, delay):
    """coincidence probability.

    Args:
        frequency1_range (array) - range of the first frequency
        dfreq1 (float) - step size of first frequency range
        frequency2_range (array) - range of the second frequency
        dfreq2 (float) - step size of second frequency range
        joint_spectral_amplitude (array((len(frequency1_range), len(frequency2_range)))) - joint spectral amplitude
        delay (array(len(delay))) - delay between photons

    Returns:
        coincidence probability array((delay))
    """
    FREQ1, FREQ2 = np.meshgrid(frequency1_range, frequency2_range)
    phase = np.exp(-1j * (FREQ1 - FREQ2) * delay)
    return 1 / 2 - 1 / 2 * np.sum(
        np.conj(joint_spectral_amplitude) * joint_spectral_amplitude.T * phase) * dfreq1 * dfreq2


general_coincidence = np.vectorize(general_coincidence_pre)


class Sellmeier:
    """Class for the index of refraction"""

    def __init__(self, A1, A2, A3, A4):
        """Constructor for the Sellmeier class.

        Args:
            (A1, A2, A3, A4) (float, float, float, float): Sellmeier parameters
        """
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.A4 = A4

    def index(self, wavelength):
        """Index of refraction.

        Args:
            wavelength (float): wavelength in micrometers
        Returns:
            array((len(wavelength)): index of refraction
        """

        return np.sqrt(self.A1 + self.A2 / (wavelength ** 2 - self.A3) - self.A4 * wavelength ** 2)


def wave_number(sellmeier, frequency):
    """wave-number calculator.

    Args:
        sellmeier (class): sellmeier class
        frequency (float): frequencies

    Returns:
        array((len(frequency)): corresponding wave-number
    """
    # calculate wavelength in units of micrometers
    wavelength = 2 * np.pi * speed_of_light * 1e6 / frequency

    # return the wave-number
    return frequency * sellmeier.index(wavelength) / speed_of_light
