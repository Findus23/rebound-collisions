import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, plot_settings

fn = filename_from_argv()
ed = ExtraData.load(fn)

plot_settings()


def get_circle(dx, dy, dz, r):
    u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
    x = r * np.cos(u) * np.sin(v) + dx
    y = r * np.sin(u) * np.sin(v) + dy
    z = r * np.cos(v) + dz
    return x, y, z


class Arrow3D(FancyArrowPatch):
    """
    https://stackoverflow.com/a/11156353/4398037
    """

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


for collision in ed.tree.get_tree().values():
    meta: CollisionMeta = collision["meta"]

    fig: Figure = plt.figure()
    ax: Axes3D = fig.gca(projection='3d')

    i = 0
    for pos, vel, r in zip(meta.collision_positions, meta.collision_velocities, meta.collision_radii):
        print(pos, vel, r)
        circle_coords = get_circle(*pos, r)
        ax.plot_wireframe(*circle_coords, color=f"C{i}")
        a = Arrow3D(*[(p, p + v / 100000) for p, v in zip(pos, vel)], mutation_scale=20,
                    lw=1, arrowstyle="-|>", color="k")
        ax.add_artist(a)
        i += 1
    ax.set_title(f"angle={meta.input.alpha:.2f}, v/v_esc={meta.input.velocity_esc:.2f}")
    xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()]).T
    diff = min(xyzlim[0] - xyzlim[1])

    print(xyzlim)
    print(diff)
    xmin, _ = ax.get_xlim3d()
    ax.set_xlim3d((xmin, xmin - diff))
    ymin, _ = ax.get_ylim3d()
    ax.set_ylim3d((ymin, ymin - diff))
    zmin, _ = ax.get_zlim3d()
    ax.set_zlim3d((zmin, zmin - diff))
    # ax.set_ylim3d(XYZlim)
    # ax.set_zlim3d(XYZlim * 3/4)
    plt.savefig("/home/lukas/tmp/3d.pdf")
    plt.show()
