import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate

from endf_parserpy import EndfParser

Threshold = 0.1


def Nuclear_trig_constants(n, eta, sign):
    """
    This function recursively evaluates the coefficients at the bottom of p. 20 in [1] for the 
    antiderivatives f_+^(n) and f_-^(n) in Appendix A.2 of [1]

    The recursion is based on the following relations:
    sin_const(n) = sign * ((n - 1) * cos_const(n - 1) + eta * sin_const(n - 1)) / ((n - 1)^2 + eta^2)
    cos_const(n) = sign * ((n - 1) * sin_const(n - 1) - eta * cos_const(n - 1)) / ((n - 1)^2 + eta^2)
    
    Parameters
    ----------
    n: degree of polynomial
    eta: eta parameter in [1]
    sign: +1 for f_+^(n), -1 for f_-^(n)
    """
    # Base case for the recursion
    if n == 0:
        return [0, 1]
    else:
        # Recursively evaluate the sin_const for n-1
        sin_const = (
            sign
            * (
                (n - 1) * Nuclear_trig_constants(n - 1, eta, sign)[0]
                + eta * Nuclear_trig_constants(n - 1, eta, sign)[1]
            )
            / ((n - 1) ** 2 + eta**2)
        )
        # Recursively evaluate the cos_const for n-1
        cos_const = (
            sign
            * (
                (n - 1) * Nuclear_trig_constants(n - 1, eta, sign)[1]
                - eta * Nuclear_trig_constants(n - 1, eta, sign)[0]
            )
            / ((n - 1) ** 2 + eta**2)
        )
        return [sin_const, cos_const]


def Trig_Integral_Eval(
    deg, cosval, sinval, y, eta, sign
):
    """
    This function evaluates the antiderivatives f_+^(deg) and f_-^(deg) at a given value y

    Parameters
    ----------
    deg: degree of polynomial
    cosval: cosine term in the antiderivative
    sinval: sine term in the antiderivative
    y: cosine of angle in centre-of-mass frame
    eta: eta parameter in [1]
    sign: +1 for f_+^(deg), -1 for f_-^(deg)
    """
    # Get constants according to recursion defined on p. 20 in [1]
    coeff = Nuclear_trig_constants(deg, eta, sign)

    # Evaluate the antiderivative at y
    return ((1 + (y * sign)) ** (deg - 1)) * (coeff[0] * sinval + coeff[1] * cosval)


def CDF_Poly_Eval(deg, eta, c, y, sign):  # evaluation of integration by parts form
    """
    This function evalutates the integral of y^deg * f_+ and y^deg * f_- in Appendix A.2 of [1] at a given value y

    Parameters
    ----------
    deg: degree of polynomial
    eta: eta parameter in [1]
    c: polar angele of the complex coefficient a_l in [1]
    y: cosine of angle in centre-of-mass frame
    sign: +1 for f_+, -1 for f_-
    """
    # Evaluate the cosine and sine terms in the antiderivative
    cosval = np.cos(eta * np.log((1 + (y * sign)) / 2) + c)
    sinval = np.sin(eta * np.log((1 + (y * sign)) / 2) + c)

    # Evaluate the sum of the antiderivative terms according to the repeated
    # integration by parts formula in Appendix A.2 of [1]
    temp = 0
    for i in range(deg + 1):
        # Trig_Integral_Eval evaluates the antiderivatives f_+^(i+1) and f_-^(i+1) at a given value y
        temp += (
            Trig_Integral_Eval(i + 1, cosval, sinval, y, eta, sign)
            * (y ** (deg - i))
            * math.factorial(deg)
            * ((-1) ** (i))
            / math.factorial(deg - i)
        )
    temp += (
        ((-1) ** (deg + 1))
        * math.factorial(deg)
        * Trig_Integral_Eval(deg + 1, cosval, sinval, 0, eta, sign)
    )
    return temp


def Legendre_Coeff_Cal(deg):  # outputs coefficients of legendre polynomial as a list
    """
    Recursively calculates the coefficients of the Legendre polynomial of degree deg.

    Based on the recurrence relation for Legendre polynomials:
    (n + 1) P_{n+1}(x) = (2n + 1) x P_n(x) - n P_{n-1}(x)

    Parameters
    ----------
    deg: degree of polynomial

    NOTE: scipy should have a function for this

    Return
    -------
    list:
        Coefficients of the Legendre polynomial of degree deg.
    """
    # Base cases for the recursion
    if deg == 0:
        return [1]
    elif deg == 1:
        return [0, 1]
    else:
        # Get coefficients for P_{deg-1}
        temp1 = Legendre_Coeff_Cal(deg - 1)
        # Shift coefficients for multiplication by x (which increases the degree by 1)
        temp1.insert(0, 0)
        # Scale coefficients according to the first term in the recurrence relation
        for i in range(len(temp1)):
            temp1[i] = (2 * (deg - 1) + 1) * temp1[i] / deg
        # Get coefficients for P_{deg-2}
        temp2 = Legendre_Coeff_Cal(deg - 2)
        temp2.append(0)
        temp2.append(0)
        # Scale coefficients according to the second term in the recurrence relation
        for i in range(len(temp2)):
            temp2[i] = (1 - deg) * temp2[i] / deg
        temp3 = []
        for i in range(len(temp2)):
            temp3.append(temp1[i] + temp2[i])
        return temp3


def CDF_Legendre_Eval(
    deg, eta, c, y, sign
): 
    """
    Evaluation of integration by parts for Legendre polynomial

    Parameters
    ----------
    deg: degree of polynomial
    eta: eta parameter on p. 8 in [1]
    c: polar angle of the complex coefficient a_l in [1]
    y: cosine of angle in centre-of-mass frame
    sign: +1 for term f_+, -1 for term f_-) in Appendix A.2 of [1]
    """
    Legendre_temp = Legendre_Coeff_Cal(deg)
    out = 0
    for i in range(deg):
        # Evaluate the integral of y^i * f_+ and y^i * f_- in Appendix A.2 of [1] at a given value y 
        # and multiply by the corresponding Legendre polynomial coefficient
        out += Legendre_temp[i] * CDF_Poly_Eval(deg, eta, c, y, sign)
    return out


def Legendre_int(deg, y):
    """
    Evaluation of Legendre polynomial integral

    Parameters
    ----------
    deg: degree of polynomial
    y: value at which to evaluate the integral of the Legendre polynomial
    """
    Legendre_temp = Legendre_Coeff_Cal(deg)
    Total = 0
    for i in range(len(Legendre_temp)):
        # antiderivative of y^i is y^(i+1)/(i+1)
        # ??? Shouldn't this be multiplied by the corresponding Legendre polynomial coefficient?
        Total += (y ** (i + 1)) / (i + 1)
    return Total


def Complex_Polar(re, im):
    """
    Convert complex number to polar form

    Parameters
    ----------
    re: real part of complex number
    im: imaginary part of complex number

    Return
    ------
    list:
        r: magnitude of complex number
        theta: angle of complex number in radians
    """
    r = np.sqrt(re**2 + im**2)
    theta = np.arctan2(im, re)
    return [r, theta]


def eta_eval(Z, z, E, m):
    """
    Compute eta on p. 142 in [4] for hydrogen
    See also p. 8 in [1] for the definition of eta
    """

    # 0.0496 = u * alpha^2 where alpha approx 1/137 is the inverse fine structure constant
    #  and u = 931.494 MeV/c^2 is the atomic mass unit
    # ? Not sure about order of magnitude. Certainly correct up to a factor 10^something.
    sqrtterm = 0.0496 * m / (E * 2)
    out = Z * z * np.sqrt(sqrtterm)
    return out


def Total_Rate_calc(
    Z, z, E, A, a, y, two_s, a_coef_re, a_coef_im, b_coef
):
    """
    Total rate between 0 and y

    Parameters
    ----------
    Z: atomic number of target particle
    z: atomic number of incident particle (proton)
    E: energy
    A: atomic mass of target particle
    a: atomic mass of incident particle (proton)
    y: cosine of angle in centre-of-mass frame
    two_s: 2 * spin of incident particle 
    a_coef_re: real part of coefficients a_l in Eq. (11) in [1]
    a_coef_im: imaginary part of coefficients a_l in Eq. (11) in [1]
    b_coef: coefficients b_l in Eq. (11) in [1]
    """
    A_ratio = A / a # ? This is just 'A' in [1] but in our case a = 1 (proton)
    eta_val = eta_eval(Z, z, E, a)

    # This is 2*eta^2/k^2 on p. 9 in [1] above eq.(11)
    # See also p. 142 in [4] for the definition of eta
    # For fundamental constants, see p. 386, Appendix H.2 Table 1 in [4]
    # I think (1e-2) * 4.16 = 2 * (h * c * alpha)^2
    # ? Not sure about order of magnitude. Certainly correct up to a factor 10^something.
    Rutherford_Coef = ((1e-2) * 4.16 * ((Z * z) * (1 + A_ratio)) ** 2) / (
        (2 * A_ratio * E) ** 2
    )
    Rate = 0

    # This is the integral of Rutherford's formula from 0 to y (see above eq.(11) on p. 9 in [1])
    # ? I tried to sanity check this by computing the derivative of this expression with respect to y and comparing it to [1].
    # ? I get an additional factor of 2 due to 2*eta_val. The derivative would agree with [1] if it was eta_val instead of 2*eta_val.
    Rate += Rutherford_Coef * (
        (y / (1 - y**2))
        + ((((-1) ** two_s) / (two_s + 1)))
        * math.sin(2 * eta_val * math.log((1 + y) / (1 - y)))
        / (2 * eta_val)
    )
    # Second sum in Eq. (11) in [1] integrated
    Nuclear_term_1 = 0
    for i in range(len(a_coef_im)):
        # Use polar form to rewrite Re(a_l * exp(i eta ...) as r * cos(theta + eta ...)
        polar_temp = Complex_Polar(a_coef_re[i], a_coef_im[i])
        int_temp = 0
        # f_+
        int_temp += ((-1) ** i) * CDF_Legendre_Eval(i, eta_val, polar_temp[1], y, 1)
        # f_-
        int_temp += CDF_Legendre_Eval(i, eta_val, polar_temp[1], y, -1)
        Nuclear_term_1 += polar_temp[0] * (2 * i + 1) * int_temp / 2
    Nuclear_term_1 *= 2 * eta_val
    Rate -= Nuclear_term_1
    Nuclear_term_2 = 0
    for i in range(len(a_coef_im)):
        # Frist sum in Eq. (11) in [1] integrated
        Nuclear_term_2 += (4 * i + 1) * b_coef[i] * Legendre_int(2 * i, y) / 2
    Rate += Nuclear_term_2
    Rate *= 2 * np.pi
    return Rate

# ? Never used (function call commented out in main code)
def Threshold_Truncator(
    ang_data, Rate_data, Threshold
):  # truncates CDF at both ends to avoid explosion in rutherford crosseciton
    i = -1
    while ang_data[i] <= Threshold:
        i -= 1
    if i == -1:
        print("error end point does not reach threshold")
    Distance = (Threshold - ang_data[(i + 1)]) / (ang_data[i] - ang_data[(i + 1)])
    New_CDF_Value = (1 - Distance) * Rate_data[(i + 1)] + Distance * Rate_data[i]
    if i < -2:
        ang_data = ang_data[: (i + 2)]
        Rate_data = Rate_data[: (i + 2)]
    ang_data[-1] = Threshold
    Rate_data[-1] = New_CDF_Value
    return ang_data, Rate_data


def ErrorCheck(Data, sign):  # checks for monotonicity in CDF and CM to Lab frame
    """
    Check for monotonicity in CDF and CM to Lab frame
    Parameters
    ----------
    Data: list of data to check
    sign: +1 for CDF (data non-decreasing with increasing index), 
          -1 for CM to Lab frame (data strictly decreasing with increasing index)
    """
    if sign == 1:
        for i in range(1, len(Data)):
            if Data[i] < Data[i - 1]:
                print("Error in CDF monotonicity")
    if sign == -1:
        for i in range(1, len(Data)):
            if Data[i] >= Data[i - 1]:
                print("Error in CM to Lab monotonicity")
                print(Data[i], Data[i - 1])
    if sign != -1 and sign != 1:
        print("sign value error")
    return 0


def CM_to_Lab_Frame(ang0, E, A, a):
    """
    Convert angle from centre-of-mass to lab frame

    Note that this is identical with the function of the same name 
    in Elastic_cross_section_generator.py.

    References:
        [1] https://doi.org/10.1088/1361-6560/ae5586


    Parameters
    ----------
    ang0: cosine of angle in centre-of-mass frame
    E: energy
    A: atomic mass of target particle
    a: atomic mass of incident particle (proton)

    Return
    ------
    Angle in lab frame
    """
    ang = math.acos(ang0)
    mp = a * 938.346  # mass of proton * c^2, MeV
    mn = A * 938.346  # mass of colliding nucleus * c^2, MeV
    p = math.sqrt((E) * (E + 2 * mp)) # p.20 in [1], note: E_1L = E + mp
    u = p / (E + mp + mn)  # ? typo should be 2mp
    g = 1 / math.sqrt(1 - u * u) # gamma_u
    e = E + mp # E_1L = E'_1L in elastic case
    v_ratio = u * (e - u * p) / (p - u * e) # u/v'_1C with v'_1c as on p. 20 in [1]

    # Eq. 13 in [1] 
    if np.fabs((g * (math.cos(ang) + v_ratio))) == 0:
        out = math.pi / 2
    else:
        out = math.atan(math.sin(ang) / (g * (math.cos(ang) + v_ratio)))
    if out < 0:
        out += math.pi
    return out


parser = EndfParser()
endf_dict = parser.parsefile("Hydrogen_elastic_data.txt")
density = {}
# Uniform spacing for cos(angle) between 0 and 0.9
angle_discretize = np.linspace(0, 0.9, 100)
angle_discretize_end = []
for i in range(
    600
):  # Setup interpolation points to match growth of rutherford crossection near boundary
    angle_discretize_end.append(np.sqrt((5.5 + 5 * i - 1) / (5.5 + 5 * i)))
angle_discretize_pos = np.concatenate((angle_discretize, angle_discretize_end))
Energy_discrete = []
Rate_discrete = []
for i in endf_dict[6][2]["subsection"][1]["E"].keys():
    print(i)
    Total_rate = []
    Real_coef = []
    Imaginary_coef = []
    b_coef = []
    Energy = endf_dict[6][2]["subsection"][1]["E"][i] * 1e-6
    if Energy >= 1:
        for r in range(8, 22):
            # Real and imaginary parts of coefficients a_l in Eq. (11) in [1]
            if r % 2 == 0:
                Real_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
            else:
                Imaginary_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
        # Coefficients b_l in Eq. (11) in [1]
        for r in range(1, 8):
            b_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
        # Compute the integral of scattering cross section from 0 to angle_discretize_pos[j] for each j
        for ang in angle_discretize_pos:
            Total_rate.append(
                Total_Rate_calc(
                    1, 1, Energy, 1, 1, ang, 1, Real_coef, Imaginary_coef, b_coef
                )
            )
        #  CM_frame_angle, Total_rate = Threshold_Truncator(CM_frame_angle, Total_rate,Threshold)
        ErrorCheck(Total_rate, 1) # Check for monotonicity in CDF
        Energy_discrete.append(Energy)
        max_rate = Total_rate[-1] # ? Reassigned before use, so this line is unnecessary
        Rate_discrete.append(max_rate) # ? Never used??
        Total_rate_True = []
        angle_discrete_True = []
        Bool_temp = True # ? Not used

        # The scattering cross section is symmetric about 0.
        # So the integral from -1+delta to y is equal to 
        # the integral from 0 to 1-delta plus sgn(y) times the integral from 0 to |y|.
        for i in range(len(Total_rate)):
            if True: # ? Not necessary
                # Get angle in lab frame for negative angles and insert into list
                angle_discrete_True.insert(
                    0, CM_to_Lab_Frame(-angle_discretize_pos[i], Energy, 1, 1)
                )
                # Insert - integral from 0 to |y| for negative angles
                Total_rate_True.insert(0, -Total_rate[i])
            
            # Get angle in lab frame for positive angles and insert into list
            angle_discrete_True.append(
                CM_to_Lab_Frame(angle_discretize_pos[i], Energy, 1, 1)
            )
            # Insert integral from 0 to |y| for positive angles
            Total_rate_True.append(Total_rate[i])

        # Retrieve integral of maximal angle noting that it was multplied by -1 
        max_rate = -Total_rate_True[0] # ? would it be easier to just use max_rate = Total_rate[-1] here?
        for i in range(len(Total_rate_True)):
            Total_rate_True[i] += max_rate # add integral of maximal angle to all values to get integral from -1+delta to y
        density[Energy] = [angle_discrete_True, Total_rate_True]

# Write data to file
# First line: energy
# Remaining lines alternate between:
# - exit angles between pi and 0
# - corresponding integrals fo scattering cross sections
f = open("hydrogen_el_ruth_cross_sec.txt", "w")
tmp = " ".join([str(z) for z in Energy_discrete])
f.write(tmp + "\n")
bool = True
for key in density.keys():
    tmp = " ".join([str(z) for z in density[key][0]])
    f.write(tmp + "\n")
    tmp = " ".join([str(z) for z in density[key][1]])
    f.write(tmp + "\n")
    bool = False
f.close()
