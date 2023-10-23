from math import pi, fmod, degrees, radians
import cmath

def sign(in_val) -> int:
    ''' Return a value based only on the numerical sign or zero of
    the input value'''
    if in_val == 0:
        return int(0)
    if in_val > 0:
        return int(1)
    if in_val < 0:
        return int(-1)

def float_to_phase(inTheta: float):
    return cmath.phase(cmath.rect(1.0,inTheta)) #(-π, π]

def normalize_angle(theta:float):
    # keep chiping away as needed
    while theta < 2*-pi or theta > 2*pi:
        theta = theta - sign(theta)*pi
    return theta

def normalize_bearing(theta:float):
    # keep chiping away as needed
    while theta < -pi or theta > pi:
        theta = theta - sign(theta)*2*pi
    return theta

def main():
    for deg in range(-405, 405, 45):
        r = radians(deg) 
        n_a = normalize_angle(r)
        n_b = normalize_bearing(r)
        sum_abs = abs(n_a) + abs(n_b) 
        print(f"{deg}\t{n_a/pi:+f} π\t{n_b/pi:+f}" )


if __name__ == '__main__':
    main()
