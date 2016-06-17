#!/usr/bin/env python

"""Plot the gravitational field on the surface of a disk"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from numpy import pi

default_R = 6000.      # radius of the disk (km)
default_H = 200.       # half-thickness of the disk (km)
default_W = 0.         # angular speed (rotation) in rad/s
earth_W = 2.*pi / (3600*24)

default_resolution = 30 # km per pixel

default_r = default_R/default_resolution
default_h = default_H/default_resolution

default_rho = 5e3      # volumetric mass (.10^3 kg/m3)
G = 6.6740831e-11     # gravitational constant


def myrhof(r, h=default_H):
    return (h**2 + r**2)**(1.5) / h


def dist(x, y, xc, yc):
    # return np.linalg.norm([x-xc, y-yc])
    return np.sqrt((x-xc)**2 + (y-yc)**2)


def calc_dist_from_center(M=None, shape=None, center=None, z=0):
    if not shape:
        cx = cy = int(M)/2
        Mx = np.array([range(M)] * M)
        My = Mx.T
    else:
        cx = int(shape[0])/2 + 1
        cy = int(shape[1])/2 + 1
        Mx = np.array([range(shape[0])] * shape[1])
        My = np.array([range(shape[1])] * shape[0])
    if center:
        cx, cy = center
    coords = [Mx - cx, My - cy]
    if z:
        coords.append(np.ones((M,M)) * z)
    distM = np.linalg.norm(coords, axis=0)
    return distM

def calc_rotation_acceleration(r, w=earth_W):
    """Given the constant angular speed w and the radius of your location r, 
    compute the acceleration"""
    return r * w**2


def fill_mass(M, rhoc=default_rho, rhof=None, resolution=default_resolution):
    """rhoc : rho is given as a constant
    rhof: rho is a function of r"""
    distM = calc_dist_from_center(M) * default_resolution
    R = distM[0, int(M)/2]
    if rhoc:
        planar = np.ones((M, M)) * rhoc
    if rhof:
        planar = rhof(distM)
    planar[distM > R] = 0
    return planar, distM


def apply_gravitation_loc(loc, planar, h=default_h,
                          resolution=default_resolution, z=None):
    if z is None:
        z = h
    M = planar.shape[0]
    R = M/2
    distloc = calc_dist_from_center(M, center=(loc, R), z=z) * resolution
    # distance in terms of x : ux directed towards the center
    # signed distances:
    distx = (np.array([range(M)]*M, dtype=float) - loc) * resolution
    disty = (np.array([range(M)]*M, dtype=float).T - R) * resolution
    distz = z * resolution

    dV = 2*h * 1.**2 * resolution**3
    dg = G * dV * planar / (distloc**2) # I am getting into trouble when dist=0
    # gravitational field
    gx = (dg * distx/distloc).sum()
    gy = (dg * disty/distloc).sum() # gy should always be zero
    gz = - (dg * distz/distloc).sum() # with z axis pointed up.
    return gx, gy, gz


def apply_gravitation(planar, step=10, resolution=default_resolution):
    Ru = planar.shape[0] / 2 # in pixel unit of the matrix 
    g = np.zeros((3, Ru))
    for loc in range(0, Ru, step):
        g[:, loc] = apply_gravitation_loc(loc, planar, resolution=resolution)
    return g

#g_side = apply_gravitation_loc(-1, planar, h=0, resolution=resolution)
def apply_rotation(g, w=earth_W, resolution=default_resolution, copy=True):
    r = np.arange(g.shape[1])[::-1] * resolution
    rot = r * (w**2)
    gr = g.copy() if copy else g
    gr[0,:] -= rot
    return gr


def draw_acceleration(g, ax, step=10, h=default_h, g_side=None,
                     resolution=default_resolution, prop=False):
    if prop:
        ax.set_aspect('equal')
    # scaling
    g_norm = np.linalg.norm(g, axis=0)
    scale_f = g_norm.max() / step
    scaled_g = g/scale_f
    max_gx = scaled_g[0].max()
    #max_gz = scaled_g[2].max()
    min_gz = scaled_g[2].min()
    xmax = g.shape[1]
    for j in range(0, xmax, step):
        gx, gy, gz = scaled_g[:, j]
        if not gx==gy==gz==0:
            ax.arrow(j-gx, 0-gz, gx, gz, length_includes_head=True,#) #,
                     head_width=step*0.1)
    if g_side is not None:
        gs_x, gs_y, gs_z = g_side
        gs_x /= scale_f
        gs_y /= scale_f
        gs_z /= scale_f
        ax.arrow(0-gs_x, -h - gs_z, gs_x, gs_z, length_includes_head=True,#) #,
                 head_width=step*0.1)
        max_gx = max(max_gx, gs_x)

    # draw the disk
    ax.add_patch(patches.Rectangle((0, 0), xmax, -2*default_h))
    ax.set_xlim(0 - 1.1*max_gx, xmax)
    ax.set_ylim(-2.1*default_h, -min_gz+step) # max_gz/scale_f
    ax.set_xticklabels(ax.get_xticks() * resolution)
    ax.set_yticklabels(ax.get_yticks() * resolution)
    ax.set_xlabel("distance from the edge to the center (km)")
    ax.set_ylabel("altitude (km)")
    ax.set_title("Gravitational field")
    # TODO: fig.add_axes on top of the patch. Imshow density.



if __name__ == '__main__':
    # Params:

    R = 6000      # radius of the disk (km)
    H = 200       # half-thickness of the disk (km)
    W = 0         # angular speed (rotation) in rad/s

    rho = 5.      # volumetric mass (.10^3 kg/m3)
    G = 6.6740831e-11 # gravitational constant

    # we represent the circle by a matrix (need to be sufficiently large
    # because of the use of "pixels")

    M = 2*R + 1
    planar = np.zeros((M, M))
    # center of the planar:
    cx = R + 1

