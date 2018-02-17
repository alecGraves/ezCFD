'''
Simple examples for fluid convection simulation in 1 dimension
'''
from numba import jit
import numpy as np


# @numba.jit(nopython=True)
@jit
def update(u, dt, dx, c, simSize):
    '''
    update the current waveform, overwriting old waveform
    Parameters
    ----------
    u : list or flat numpy array of floats
        The waveform to be updated (overwritten)
    dx : float
        Distance between x (spatial) values
    dt : float
        Distance between simulation time steps
    c : float
        The wavespeed (default 1.0)
    '''
    un = u.copy()
    scale = c*dt/dx
    for i in range(1, simSize):
        un[i] = u[i] - scale*(u[i]-u[i-1])
    return un



# @numba.jit(nopython=True)
@jit
def updateInviscid(u, dt, dx, simSize):
    '''
    update the current waveform, overwriting old waveform
    uses Inviscid Burgers
        -> Can gernerate discontinuous solutions from
            smooth initial conditions (similar to shocks 
            in supersonic flows).
    Parameters
    ----------
    u : list or flat numpy array of floats
        The waveform to be updated (overwritten)
    dx : float
        Distance between x (spatial) values
    dt : float
        Distance between simulation time steps
    '''
    un = u.copy()
    scale = dt/dx
    # scale = 0.5
    for i in range(1, simSize):
        un[i] = u[i] - scale*u[i]*(u[i]-u[i-1])
    return un

class LinearConvectionSim(object):
    '''
    Class to compute linear convection using 'du/dt + c(du/dx) = 0'
    Solution uses forward difference for numerical (du/dt)
        and backward difference for numerical (du/dx).
    numerically solves for next state in terms of previous state.
    
    
    '''
    def __init__(self, u, x, t, c=1.0):
        '''
        pass in initial wave, x domain and t domain
        Parameters
        ----------
        u : list or flat numpy array of floats
            The initial waveform for the simulation.
        x : list or flat numpy array of floats
            The x indices for every u, evenly spaced (linspace).
        t : list or flat numpy array of floats
            The timesteps to use for the similation (evenly spaced).
        c : float
            The wavespeed (default 1.0)
        '''
        # given:
        self.simSize = len(u)
        self.u = u
        self.x = x
        self.dx = x[1] - x[0]
        self.t = t
        self.dt = t[1] - t[0]
        self.c = c
    
    def step(self, mode='inviscid'):
        '''
        generator which yields updates to the simulation
        Parameters
        ----------
        mode : str
            The simulation mode, either "inviscid" or "regular"
        '''
        if 'inviscid' in mode:
            self.u = updateInviscid(self.u, self.dt, self.dx, self.simSize)
        else:
            self.u = update(self.u, self.dt, self.dx, self.c, self.simSize)

def _visualize_test():
    import matplotlib.pyplot as plt
    import time

    x = np.linspace(0, 2, 41)
    t = np.linspace(0, 0.5, 25)
    u = np.ones((41,))
    u[(x>0.5) & (x<1)] = 2.

    l = LinearConvectionSim(u, x, t)

    plt.ion()
    hl, = plt.plot(x, l.u)
    plt.show()
    plt.pause(1e-17)

    for _time in t:
        l.step('inviscid')
        hl.set_ydata(l.u)
        plt.title(str(int(_time*1000)) + ' ms')
        plt.draw()
        plt.pause(1e-17)
        time.sleep(0.025)
    
    l = LinearConvectionSim(u, x, t)

    for _time in t:
        l.step('regular')
        hl.set_ydata(l.u)
        plt.title(str(int(_time*1000)) + ' ms')
        plt.draw()
        plt.pause(1e-17)
        time.sleep(0.025)


if __name__ == "__main__":
    _visualize_test()


