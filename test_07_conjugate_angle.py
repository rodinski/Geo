import Geo


for i in range(-360, 360, 10):
    a = Geo.Angle(i, unit='deg')
    print(f"{a=}    {Geo.conjugate_angle(a)=}")
