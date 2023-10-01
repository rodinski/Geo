import Geo
import matplotlib.pyplot as plt
from matplotlib import patches

pts = []
seg = []
cur = []
ray = []
pts.append(Geo.Point(0, 0))
pts.append(Geo.Point(88, 22))
seg.append(Geo.Segment(pts[0], pts[1]))

seg.append(seg[0].copy_parallel(-110))
seg.append(seg[0].copy_parallel(+110))
#print(seg[-1])
print(len(seg))

cur.append(Geo.Curve(pts[0], pts[1], Geo.Angle(5, unit='deg')))
cur.append(cur[0].copy_parallel(-30))
cur.append(cur[0].copy_parallel(-60))
cur.append(cur[0].copy_parallel(-80))

ray.append(Geo.Ray(Geo.Point(12,3), 0.14))
ray.append(ray[0].copy_parallel(-10))
ray.append(ray[0].copy_parallel(1))


fig, ax = plt.subplots()
for i in seg:
    print(i)
    ax.add_patch(i.patch())
for curve in cur:
    for p in curve.patch_all():
        ax.add_patch(p)
for r in ray:
    print(r.patch())
    ax.add_patch(r.patch())
    ax.scatter(r.Point.X, r.Point.Y)

plt.axis('scaled')
plt.show()


