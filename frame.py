"""
Module frame renders images of the 2D system, with particles rendered as
circles and arrows indicating their displacement.

(modified from
https://github.com/yketa/active_particles/tree/master/analysis/frame.py)

Environment modes
-----------------
MODE : string
    Plotting mode.
     _________________________________________________________________________
    | Mode           | Arrow direction        | Arrow length | Particle color |
    |________________|________________________|______________|________________|
    | 'orientation'  | Self-propulsion        | Relative to  | Orientation    |
    |                | direction              | diameter     |                |
    |________________|________________________|______________|________________|
    DEFAULT: orientation
PLOT : bool
    Plot single frame.
    DEFAULT: False
MOVIE : bool
    Make movie out of several plotted frames.
    DEFAULT: False
SHOW [PLOT mode] : bool
    Show figure.
    DEFAULT: False
SAVE [PLOT mode] : bool
    Save figure.
    DEFAULT: False
SUPTITLE : bool
    Display suptitle on figures.
    DEFAULT: True

Environment parameters
----------------------
DAT_FILE : string
    Data file.
    DEFAULT: CURRENT_DIRECTORY/out.dat
INITIAL_FRAME : int
    [PLOT mode] Frame to render.
    [MOVIE mode] Initial frame to render.
    NOTE: FRAME < 0 will be interpreted as the frame to render being the middle
          frame of the simulation.
    DEFAULT: -1
FINAL_FRAME [MOVIE mode] : int
    Final movie frame.
    DEFAULT: Final simulation frame.
FRAME_PERIOD [MOVIE mode] : int
    Frame rendering period.
    DEFAULT: active_work._frame_per
FRAME_MAXIMUM : int
    Maximum number of frames.
    DEFAULT: active_work.frame._frame_max
DT : int
    Lag time for displacement.
    NOTE: We consider time step in fixed in simulations.
    NOTE: [PLOT mode] DT < 0 will be interpreted as a lag time corresponding to
                      the total number of simulation frames - FRAME + TIME.
          [MOVIE mode] DT < 0 will be interpreted as a lag time corresponding
                       to the minimum distance between frames.
    DEFAULT: -1
BOX_SIZE : float
    Length of the square box to render.
    DEFAULT: simulation box size
X_ZERO : float
    1st coordinate of the centre of the square box to render.
    DEFAULT: 0
Y_ZERO : float
    2nd coordinate of the centre of the square box to render.
    DEFAULT: 0
V_MIN : float
    Minimum value of the colorbar.
    DEFAULT: 10^{E(log(||\\vec{u}||))-2*V(log(||\\vec{u}||))}
    NOTE: Colorbar is represented in logarithmic scale so V_MIN > 0.
V_MAX : float
    Maximum value of the colorbar.
    DEFAULT: 10^{E(log(||\\vec{u}||))+2*V(log(||\\vec{u}||))}
    NOTE: Colorbar is represented in logarithmic scale so V_MAX > 0.
LABEL : bool
    Write indexes of particles in circles.
    DEFAULT: False
ARROW_WIDTH : float
    Width of the arrows.
    DEFAULT: active_work.frame._arrow_width
HEAD_WIDTH : float
    Width of the arrows' head.
    DEFAULT: active_work.frame._arrow_head_width
HEAD_LENGTH : float
    Length of the arrows' head.
    DEFAULT: active_work.frame._arrow_head_length
FRAME_VERTICAL_SIZE : float
    Vertical size of the frame (in inches).
    DEFAULT: active_work.frame._frame_ver
FRAME_HORIZONTAL_SIZE : float
    Horizontal size of the frame (in inches).
    DEFAULT: active_work.frame._frame_hor
FRAME_DEFINITION [SAVE mode] : float
    Definition of image (in dots per inches (dpi)).
    DEFAULT: active_work.frame._frame_def
FONT_SIZE : int
    Font size.
    DEFAULT: active_work.frame._font_size
PAD : float
    Separation between label and colormap.
    DEFAULT: active_work.frame._colormap_label_pad
FIGURE_NAME [SAVE or MOVIE mode] : string
    Custom figure name.
    DEFAULT: out.eps [figure] out.mp4 [movie]
MOVIE_DIR [MOVIE mode] : string
    Custom movie directory.
    DEFAULT: out.movie

Output
------
> Prints execution time.
[PLOT mode]
> Plots system according to plotting mode.
[MOVIE mode]
> Creates movie from frame concatenation according to plotting mode.
[SHOW mode]
> Displays figure.
[SAVE mode]
> Saves figure in current directory.
"""

from active_work.init import get_env, mkdir
from active_work.read import Dat
from active_work.maths import normalise1D, amplogwidth

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

import sys

from math import ceil

import numpy as np
np.seterr(divide='ignore')

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as ColorsNormalise
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

from datetime import datetime

from collections import OrderedDict

import subprocess

# DEFAULT VARIABLES

_frame_per = 1      # default frame rendering period
_frame_max = 1000   # default maximum number of frames

_frame_ver = 12 # default vertical size of the frames (in inches)
_frame_hor = 16 # default horizontal size of the frames (in inches)
_frame_def = 80 # default definition of images (in dots per inches (dpi))

_arrow_width = 1e-3                         # default width of the arrows
_arrow_head_width = _arrow_width*3e2        # default width of the arrows' head
_arrow_head_length = _arrow_head_width*1.5  # default length of the arrows' head

_font_size = 15 # font size

_colormap_label_pad = 30    # default separation between label and colormap

# FUNCTIONS AND CLASSES

class _Frame:
    """
    This class is designed as the superclass of all other plotting classes
    specific to each mode. It initialises the figure and provides methods to
    plot circles representing particles and arrows at the particles' positions,
    and to add a colorbar.
    """

    def __init__(self, dat, frame, box_size, centre,
        arrow_width=_arrow_width,
        arrow_head_width=_arrow_head_width,
        arrow_head_length=_arrow_head_length):
        """
        Initialises figure.

        Parameters
        ----------
        dat : active_particles.read.Dat
    		Data object.
        frame : int
            Frame to render.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        """

        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim([-1.1*box_size/2, 1.1*box_size/2])
        self.ax.set_xlabel(r'$x$')
        self.ax.set_ylim([-1.1*box_size/2, 1.1*box_size/2])
        self.ax.set_ylabel(r'$y$')
        self.ax.set_aspect('equal')

        self.positions = dat.getPositions(frame, centre=centre) # particles' positions at frame frame with centre as centre of frame
        self.diameters = np.full((dat.N,), fill_value=1)        # particles' diameters at frame frame

        self.particles = [particle for particle in range(len(self.positions))
            if (np.abs(self.positions[particle]) <= box_size/2).all()]  # particles inside box of centre centre and length box_size

        self.arrow_width = arrow_width
        self.arrow_head_width = arrow_head_width
        self.arrow_head_length = arrow_head_length

    def __del__(self):
        """
        Closes figure.
        """

        plt.close(self.fig)

    def draw_circle(self, particle, color='black', fill=False, label=False):
        """
        Draws circle at particle's position with particle's diameter.

        Parameters
        ----------
        particle : int
            Particle index.
        color : any matplotlib color
            Circle color. (default: 'black')
        fill : bool
            Filling the circle with same color. (default: False)
        label : bool
            Write indexes of particles in circles. (default: False)
        """

        circle = plt.Circle(self.positions[particle],
            self.diameters[particle]/2, color=color, fill=fill, zorder=0)   # circle artist representing particle
        self.ax.add_artist(circle)
        if label:
            self.ax.annotate(
                "%i" % particle, xy=self.positions[particle], ha="center")

    def draw_arrow(self, particle, dx, dy, color='black'):
        """
        Draws arrow starting from particle's position.

        Parameters
        ----------
        particle : int
            Particle index.
        dx : float
            Arrow length in x-direction.
        dy : float
            Arrow length in y-direction.
        color : any matplotlib color
            Arrow color. (default: 'black')
        """

        length = np.sqrt(dx**2 + dy**2) # length of arrow
        if length == 0: return
        self.ax.arrow(*self.positions[particle], dx, dy, color=color,
            width=length*self.arrow_width,
            head_width=length*self.arrow_head_width,
            head_length=length*self.arrow_head_length, zorder=1)

    def colorbar(self, vmin, vmax, cmap=plt.cm.jet):
        """
        Adds colorbar to plot.

        Parameters
        ----------
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        cmap : matplotlib colorbar
            Matplotlib colorbar to be used. (default: matplotlib.pyplot.cm.jet)
        """

        vNorm = ColorsNormalise(vmin=vmin, vmax=vmax)
        self.scalarMap = ScalarMappable(norm=vNorm, cmap=cmap)

        self.colormap_ax = make_axes_locatable(self.ax).append_axes('right',
            size='5%', pad=0.05)
        self.colormap = mpl.colorbar.ColorbarBase(self.colormap_ax, cmap=cmap,
            norm=vNorm, orientation='vertical')

class Orientation(_Frame):
    """
    Plotting class specific to 'orientation' mode.
    """

    def __init__(self, dat, frame, box_size, centre,
        arrow_width=_arrow_width,
        arrow_head_width=_arrow_head_width,
        arrow_head_length=_arrow_head_length,
        pad=_colormap_label_pad, dt=0,
        label=False, **kwargs):
        """
        Initialises and plots figure.

        Parameters
        ----------
        dat : active_particles.read.Dat
    		Data object.
        frame : int
            Frame to render.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        pad : float
            Separation between label and colormap.
            (default: active_work.frame._colormap_label_pad)
        dt : int
            Lag time for displacement. (default=0)
        label : bool
            Write indexes of particles in circles. (default: False)
        """

        super().__init__(dat, frame, box_size, centre,
            arrow_width=arrow_width,
            arrow_head_width=arrow_head_width,
            arrow_head_length=arrow_head_length)    # initialise superclass

        self.orientations = (
            dat.getOrientations(frame, *self.particles)%(2*np.pi))  # particles' orientations at frames

        self.colorbar(0, 2, cmap=plt.cm.hsv)                                    # add colorbar to figure
        self.colormap.set_label(r'$\theta_i/\pi$', labelpad=pad, rotation=270)  # colorbar legend

        self.label = label  # write labels

        self.draw()

    def draw_arrow(self, particle, orientation, color='black'):
        """
        Draws arrow in particle.

        Parameters
        ----------
        particle : int
            Particle index.
        orientation : float
            Orientation of particle in radians.
        color : any matplotlib color
            Arrow color. (default: 'black')
        """

        direction = np.array([np.cos(orientation), np.sin(orientation)])
        length = self.diameters[particle]*0.75
        self.ax.arrow(
            *(self.positions[particle] - direction*length/(2*np.sqrt(2))),
            *direction*length/np.sqrt(2),
            color=color,
            width=length*self.arrow_width,
            head_width=length*self.arrow_head_width,
            head_length=length*self.arrow_head_length, zorder=1,
            length_includes_head=True)

    def draw(self):
        """
        Plots figure.
        """

        for particle, orientation in zip(self.particles,
            self.orientations):                     # for particle and particle's displacement in rendered box
            self.draw_circle(particle,
                color=self.scalarMap.to_rgba(orientation/np.pi), fill=True,
                label=self.label)                   # draw particle circle with color corresponding to displacement amplitude
            self.draw_arrow(particle, orientation)  # draw displacement direction arrow

# SCRIPT

if __name__ == '__main__':  # executing as script

    startTime = datetime.now()

    # VARIABLE DEFINITIONS

    mode = get_env('MODE', default='orientation')           # plotting mode
    if mode == 'orientation':
        plotting_object = Orientation
    else: raise ValueError('Mode %s is not known.' % mode)  # mode is not known

    dat_file = get_env('DAT_FILE', default=joinpath(getcwd(), 'out.dat'))   # data file
    dat = Dat(dat_file)                                                     # data object

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)  # initial frame to render

    box_size = get_env('BOX_SIZE', default=dat.L, vartype=float)    # size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
        get_env('Y_ZERO', default=0, vartype=float))                # centre of the box

    Nentries = dat.frames - 1
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame    # initial frame to draw

    # FIGURE PARAMETERS

    vmin = get_env('V_MIN', vartype=float) # minimum value of the colorbar
    vmax = get_env('V_MAX', vartype=float) # maximum value of the colorbar

    frame_hor = get_env('FRAME_HORIZONTAL_SIZE', default=_frame_hor,
        vartype=float)  # horizontal size of the frame (in inches)
    frame_ver = get_env('FRAME_VERTICAL_SIZE', default=_frame_ver,
        vartype=float)  # vertical size of the frame (in inches)
    mpl.rcParams['figure.figsize'] = (frame_hor, frame_ver)

    frame_def = get_env('FRAME_DEFINITION', default=_frame_def,
        vartype=float)                                                  # definition of image (in dots per inches (dpi))
    font_size = get_env('FONT_SIZE', default=_font_size, vartype=float) # font size
    mpl.rcParams.update({'savefig.dpi': frame_def, 'font.size': font_size})

    arrow_width = get_env('ARROW_WIDTH', default=_arrow_width,
        vartype=float)  # width of the arrows
    arrow_head_width = get_env('HEAD_WIDTH', default=_arrow_head_width,
        vartype=float)  # width of the arrows' head
    arrow_head_length = get_env('HEAD_LENGTH', default=_arrow_head_length,
        vartype=float)  # length of the arrows' head

    pad = get_env('PAD', default=_colormap_label_pad, vartype=float)    # separation between label and colormap

    # LEGEND SUPTITLE

    display_suptitle = get_env('SUPTITLE', default=True, vartype=bool)  # display suptitle

    def suptitle(frame, lag_time=None):
        """
        Returns figure suptitle.

        NOTE: Returns empty string if display_suptitle=False.

        Parameters
        ----------
        frame : int
            Index of rendered frame.
        lag_time : int
            Lag time between frames.

        Returns
        -------
        suptitle : string
            Suptitle.
        """

        if not(display_suptitle): return ''

        suptitle = (
            str(r'$N=%.2e, \phi=%1.2f, l_p/\sigma=%.2e$'
    		% (dat.N, dat.phi, dat.lp)))
        suptitle += str(r'$, L=%.3e$' % dat.L)
        if 'BOX_SIZE' in envvar:
            suptitle += str(r'$, L_{new}=%.3e$' % box_size)
        suptitle += '\n'
        if 'X_ZERO' in envvar or 'Y_ZERO' in envvar:
            suptitle += str(r'$x_0 = %.3e, y_0 = %.3e$' % centre) + '\n'
        suptitle += str(r'$t/(l_p/\sigma) = %.5e$'
            % (frame*dat.dt/dat.lp))
        if lag_time != None:
            suptitle += str(r'$, \Delta t/(l_p/\sigma) = %.5e$'
                % (lag_time*dat.dt/dat.lp))

        return suptitle

    # MODE SELECTION

    if get_env('PLOT', default=False, vartype=bool):    # PLOT mode

        Nframes = Nentries - init_frame  # number of frames available for the calculation

        figure = plotting_object(dat, init_frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length,
            pad=pad, dt=dt, vmin=vmin, vmax=vmax,
            label=get_env('LABEL', default=False, vartype=bool))
        figure.fig.suptitle(suptitle(init_frame))

        if get_env('SAVE', default=False, vartype=bool):    # SAVE mode
            figure_name, = naming_standard.filename(**attributes)
            figure.fig.savefig(joinpath(data_dir,
                get_env('FIGURE_NAME', default='out.eps')))

    if get_env('MOVIE', default=False, vartype=bool):   # MOVIE mode

        frame_fin = get_env('FINAL_FRAME', default=Nentries, vartype=int)       # final movie frame
        frame_per = get_env('FRAME_PERIOD', default=_frame_per, vartype=int)    # frame rendering period
        frame_max = get_env('FRAME_MAXIMUM', default=_frame_max, vartype=int)   # maximum number of frames

        movie_dir = get_env('MOVIE_DIR', default='out.movie')   # movie directory name
        mkdir(movie_dir)                                        # create movie directory
        mkdir(joinpath(movie_dir, 'frames'), replace=True)      # create frames directory (or replaces it if existing)

        frames = [init_frame + i*frame_per for i in range(frame_max)
            if init_frame + i*frame_per <= frame_fin]   # rendered frames

        for frame in frames:    # for rendered frames
            sys.stdout.write(
                'Frame: %d' % (frames.index(frame) + 1)
                + "/%d \r" % len(frames))

            figure = plotting_object(dat, frame, box_size, centre,
                arrow_width, arrow_head_width, arrow_head_length,
                pad=pad, dt=frame_per, vmin=vmin, vmax=vmax,
                label=get_env('LABEL', default=False, vartype=bool))    # plot frame
            figure.fig.suptitle(suptitle(frame, frame_per))

            figure.fig.savefig(joinpath(movie_dir, 'frames',
                '%010d' % frames.index(frame) + '.png'))    # save frame
            del figure                                      # delete (close) figure

        subprocess.call([
            'ffmpeg', '-r', '5', '-f', 'image2', '-s', '1280x960', '-i',
            joinpath(movie_dir , 'frames', '%10d.png'),
            '-pix_fmt', 'yuv420p', '-y',
            joinpath(movie_dir, get_env('FIGURE_NAME', default='out.mp4'))
            ])  # generate movie from frames

    # EXECUTION TIME
    print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('SHOW', default=False, vartype=bool):    # SHOW mode
        plt.show()
