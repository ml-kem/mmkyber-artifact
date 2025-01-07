load("https://bitbucket.org/malb/lwe-estimator/raw/HEAD/estimator.py")

def findMLWEdelta(bar_nu, n, d, logq):
    n = n * d
    q = 2^logq
    stddev = float(sqrt((bar_nu^2 - 1)/12))
    alpha = alphaf(sigmaf(stddev),q)
    set_verbose(1)
    L = estimate_lwe(n, alpha, q, reduction_cost_model=BKZ.enum)
    delta_enum1 = L['usvp']['delta_0']
    delta_enum2 = L['dec']['delta_0']
    delta_enum3 = L['dual']['delta_0']
    L = estimate_lwe(n, alpha, q, reduction_cost_model=BKZ.sieve)
    delta_sieve1 = L['usvp']['delta_0']
    delta_sieve2 = L['dec']['delta_0']
    delta_sieve3 = L['dual']['delta_0']
    return max(delta_enum1,delta_enum2,delta_enum3,delta_sieve1,delta_sieve2,delta_sieve3)


def find_MLWE_n(nu, d, logq, sec):
    mlwe =  1                               # dimension of the Module-LWE problem
    mlwe_hardness = 2
    # 128-bit security: 1.0045
    # 192-bit security: 1.0029
    # 256-bit security: 1.0023
    if sec == 128:
        while mlwe_hardness > 1.0045:           # increasing the mlwe dimension until MLWE provides ~ 128-bit security
            mlwe += 1
            mlwe_hardness = findMLWEdelta(nu,mlwe,d, logq)
    if sec == 192:
        while mlwe_hardness > 1.0029:           # increasing the mlwe dimension until MLWE provides ~ 192-bit security
            mlwe += 1
            mlwe_hardness = findMLWEdelta(nu,mlwe,d, logq)
    if sec == 256:
        while mlwe_hardness > 1.0023:           # increasing the mlwe dimension until MLWE provides ~ 192-bit security
            mlwe += 1
            mlwe_hardness = findMLWEdelta(nu,mlwe,d, logq)
    return mlwe

def find_sigma(d, m ,n ):
   return sqrt(2) * sqrt( ln(2 * (d * (m + n)) * (1 + 2^(1+128))) / pi)

d = 256
nu = 1
bar_nu = 3
sec = 128
logq = 24
ell = 2^10
d_v = 2
tau = sqrt(-log(2^(-129))/pi)

flag = 0
while flag == 0:
    n = find_MLWE_n(bar_nu, d, logq, sec)

    m = n
    sigma = find_sigma(d, m, n)
    sigma_1 = 2 * sqrt(sigma^2 + 1)

    B = ell * (n+m) * (d * nu)^2

    sigma_2 = 2 * sqrt(B) * sqrt(sigma^2 + 1)
    d_u = logq
    dec_error = (n+m) * d * nu * sigma_1 * tau + sigma_2 * tau + round(2^logq / 2^(d_v + 1))
    if dec_error < 2^logq / 4:
        flag = 1
        print("sec level = ", sec)
        print("mmPKE: N = ", ell)
        print("mmPKE: logq = ", logq)
        print("mmPKE: n = ", n)
        print("mmPKE: public parameter |pp| size (KB) = ", round((float(d*logq*n^2)/8/1024),3))
        print("mmPKE: public key |pk| size (KB) = ", round((float(d*logq*n)/8/1024),3))

        d_u = 0
        dec_error_1 = (n+m) * d * nu * sigma_1 * tau + sigma_2 * tau + round(2^logq / 2^(d_v + 1)) + nu * d * n * round(2^logq / 2^(d_u + 1))
        while dec_error_1 >= 2^logq / 4 and d_u <= logq:
            d_u = d_u + 1
            dec_error_1 = (n+m) * d * nu * sigma_1 * tau + sigma_2 * tau + round(2^logq / 2^(d_v + 1)) + nu * d * n * round(2^logq / 2^(d_u + 1))
        print("d_u = ", d_u)
        print("mmPKE: log (dec error') = ", round(log(dec_error_1, 2), 3))
        print("mmPKE: |ct_0| (KB) = ", round((float(d*d_u*(n))/8/1024), 3))
        print("mmPKE: |ct_i| (KB) = ", round((float(d*d_v)/8/1024), 3))
        print("mmPKE: total |ct| (KB) = ", round((float(d*d_v*ell + d*d_u*(n))/8/1024), 3))
        print("mmPKE-comp: compress rate k = ", round(float(ell * 768*8*sec/128)/float(d*d_u*n + d*d_v*ell),3) )


    else:
        logq = logq + 1
        print("logq = ", logq)


find_prime = 0
q = 2^logq
q = q + 1
while find_prime == 0:
    if is_prime(q) == True:
        find_prime = 1
        print("prime q = ", q)
        print("[log q] = ", round(log(q, 2),3))
        if dec_error >= 2^logq / 4:
            print("!!!!! q fail")
    q = q - 2^8


print("Gaussian width sigma_0 = ", round(sigma_1, 3))
print("Gaussian width sigma_1 = ", round(sigma_2, 3))
print("Standard deviation s_0 = ", round(sigma_1/sqrt(2*pi), 3))
print("Standard deviation s_1 = ", round(sigma_2/sqrt(2*pi), 3))

