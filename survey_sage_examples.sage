## Code for examples in the paper
## ``Survey: Recovering cryptographic keys from partial information, by example''
## by Gabrielle De Micheli and Nadia Heninger


# Section 4.2.1 : Lattice attacks on low-exponent RSA with bad padding

def low_exp_rsa():
    N = 0x98664cf0c9f8bbe76791440d
    a = 0x01FFFFFFFFFFFFFFFF0000 
    c = 0xeb9a3955a7b18d27adbf3a1

    R.<x> = ZZ[]

    f = (a + x)^3 - c
    print("polynomial is", expand(f) % N)
    # f_expand = x^3 + 0x5fffffffffffffffd0000*x^2 + 0x6f1c485f406ba1c069460efe*x + 0x203211880cdc43afe1c5c5f9

    f2 = f[2] % N  #0x5fffffffffffffffd0000
    f1 = f[1] % N # 0x6f1c485f406ba1c069460efe
    f0 = f[0] % N  # 0x203211880cdc43afe1c5c5f9
    
    M = matrix(4)
    X = 2^16
    M[0] = [X^3, X^2*f2, X*f1, f0]
    M[1,1] = N*X^2
    M[2,2] = N*X
    M[3,3] = N

    A = M.LLL()
    
    g = A[0][0]/X^3*x^3 + A[0][1]/X^2 *x^2 + A[0][2]/X *x + A[0][3]

    return hex(g.roots()[0][0])


# Section 4.2.2 : Factorization from consecutive bits of p

def facto_consec():

    N = 0x4d14933399708b4a5276373cb5b756f312f023c43d60b323ba24cee670f5
    a = 0x68323401cb3a10959e7bfdc0000000

    X = 2^30

    M = matrix(3)
    M[0] = [X^2, X*a, 0]
    M[1] = [0, X, a]
    M[2,2] = N 

    A = M.LLL()

    R.<x> = ZZ[]
    f = A[0][2] + A[0][1]/X *x + A[0][0]/X^2*x^2
    print("polynomial f is", f)

    r = f.roots()[0][0]

    print("f has one positive integer root", hex(r))

    return N/gcd(int(a) + r, N)


# ------------Section 4.2.4 : RSA Key recovery from middle bits of p---------------------------

def bivariate_coppersmith():
    R.<x,y> = ZZ[]
    #p = random_prime(2^164,2^163)
    #q = random_prime(2^164,2^163)
    #N = p*q
    #a = lift(mod(p,2^148)) - lift(mod(p,2^16)) 

    N = 0x3ab05d0c0694c6bd8ee9683d15039e2f738558225d7d37f4a601bcb929ccfa564804925679e2f3542b
    a = 0xc48c998771f7ca68c9788ec4bff9b40b80000
    
    X = 2^16
    Y = 2^16

    f = x + a + y*2^148
    monomial_list = (f^3).monomials()
    function_list = [f^3,f^2*y,f*(y)^2, (y)^3*N, f^2, f*(y), (y)^2*N, f, (y)*N, N]

    M = matrix(10)
    for i in range(10):
        M[i] = [R(function_list[i])(x*X,y*Y).monomial_coefficient(m) for m in monomial_list]

    scaled_monomials = [m(x/X,y/Y) for m in monomial_list]
        
    def getf(M,i):
        return sum(b*m for b,m in zip(M[i],scaled_monomials))

    A = M.LLL()

    print("N =",hex(N))
    print("a =",hex(a))

    I = ideal(getf(A,0),getf(A,1))
    print(I.groebner_basis())

    I = ideal(getf(A,i) for i in range(3))
    print(I.groebner_basis())
    
    I = ideal(getf(A,i) for i in range(9))
    print(I.groebner_basis())

   # print(getf(A,0))

   # print(getf(A,1))
    
def print_matrix():
    RR.<x,y,T,a,R,N> = ZZ[]

    f = x + T * y + a
    monomial_list = ((x+y+1)^3).monomials()
    function_list = [f^3,f^2*y,f*(y)^2, (y)^3*N, f^2, f*(y), (y)^2*N, f, (y)*N, N]

    for i in range(10):
        for m in monomial_list:
            s = str(RR(function_list[i])(x=x*R,y=y*R).coefficient(m)(x=0,y=0)).replace("*","")
            print(s,end=" & ")
        print("\\\\")
            
proof.arithmetic(False)




# Section 4.2.7 : Partial recovery of RSA dp

def rsa_dp():

    N = 0x4d14933399708b4a5276373cb5b756f312f023c43d60b323ba24cee670f5
    a = 0x25822d06984a06be5596fcc0000000
    e = 65537
    kp = 23592

    A = a + inverse_mod(e, N)*(kp-1)

    print("A is", hex(A))

    A = 0x8ffe9143aa4c189787058057a0784576848f3f28d79a83169f72a0550699112

    X = 2^30

    M = matrix(3)
    M[0] = [X^2, X*A, 0]
    M[1] = [0, X, A]
    M[2,2] = N

    B = M.LLL()

    R.<x> = ZZ[]
    f = B[0][2] + B[0][1]/X *x + B[0][0]/X^2*x^2
    print("polynomial f is", f)

    r = f.roots()[0][0]

    print("f has integer root", hex(r))

    return N/gcd(int(A) + r, N)



# Section 5.2.1 : Lattice attacks

def lattice_attacks():

    p,F,C,n,G,x = ecdsa_params()

    print("n", hex(n))
    print("p", hex(p))
    print("G", G)

    h1 = 0xae0f1d8cd0fd6dd1
    h2 = 0x8927e246fe4f3941

    sig1 = '6393e79fbfb40c9c 621ee64e65d1e938'
    r1,s1 = [Integer(f,16) for f in sig1.split()]
    sig2 = '3ea8720afa6d03c2 16fc6aa65bf241ea'
    r2,s2 = [Integer(f,16) for f in sig2.split()]

    print(hex(r2))

    t = Integer(-inverse_mod(s1,n)*s2*r1*inverse_mod(r2,n))
    u = Integer(inverse_mod(s1,n)*r1*h2*inverse_mod(r2, n)-inverse_mod(s1,n)*h1)

    K = 2^32

    M = matrix(3)
    M[0, 0] = n
    M[1] = [t, 1, 0]
    M[2] = [u, 0, K]


    A = M.BKZ()

    print(A)

    k1 = A[2][0]
    k2 = A[2][1]

    print((k1*G)[0] == r1)

    return "vector (k1, k2)", (hex(k1), hex(k2))

# ---------------  Section 5.2.3 : (EC)DSA key recovery from middle bits of the nonce k -----------------------------

def ecdsa_middle_bits():
    p,F,C,n,G,x = ecdsa_params()
    
    h1 = 0x608932fcfaa7785d
    h2 = 0xe5f8eca48ac2a45c

    k1 = 0x734450e2fd5da41c
    sig1 = '1a4adeb76b4a90e0 eba129bb2f97f7cd'
    r1,s1 = [Integer(f,16) for f in sig1.split()]
    k2 = 0x4de972930ab4a534
    sig2 = 'c4e5bec792193b51 0202d6eecb712ae3'
    r2,s2 = [Integer(f,16) for f in sig2.split()]

    a1 = lift(mod(k1,2^(64-15)))-lift(mod(k1,2^15))
    a2 = lift(mod(k2,2^(64-15)))-lift(mod(k2,2^15))

    print("a1=",hex(a1))
    print("a2=",hex(a2))
    
    b1 = lift(mod(k1,2^15))
    b2 = lift(mod(k2,2^15))
    
    c1 = 2^(-64+15)*(k1 - lift(mod(k1,2^(64-15))))
    c2 = 2^(-64+15)*(k2 - lift(mod(k2,2^(64-15))))

    t = Integer(r1*inverse_mod(s1,n)*inverse_mod(r2,n)*s2)
    u = Integer(-inverse_mod(s1,n)*h1+r1*inverse_mod(s1,n)*inverse_mod(r2,n)*h2)

    print(mod(b1+c1*2^(64-15)-t*b2-t*c2*2^(64-15)+a1-t*a2+u,n))

    M = matrix(5)
    X = 2^15
    M[0] = [X, X*2^(64-15), -X*t, -X*t*2^(64-15), a1-t*a2+u]
    M[1,1] = n*X
    M[2,2] = n*X
    M[3,3] = n*X
    M[4,4] = n

    A = M.LLL()
    
    R.<x1,y1,x2,y2> = ZZ[]
    
    def getf(M,i):
        return M[i,0]/X*x1+M[i,1]/X*y1+M[i,2]/X*x2+M[i,3]/X*y2+M[i,4]

    I = ideal(getf(A,i) for i in range(4))
    return I.groebner_basis()
    

# ---------------  Section 6.2 : Most significant bits of finite field Diffie-Hellman shared secret --------------------

def dh_msb():
    p,g = dh_params()

    c = 0x12d0dca5769537c3cd47d8f9042f7497
    d = 0x45fb9bfbdcbead5616aacc7b0f879ae4

    DH1 = lift(mod(g,p)^(d*c))

    r = 0x56e112dac14f4a4cc02951414aa43a38

    DH2 = lift(mod(g,p)^((d+r)*c))

    a1 = DH1 - lift(mod(DH1,2^63))
    a2 = DH2 - lift(mod(DH2,2^63))

    b1 = DH1-a1
    b2 = DH2-a2
    
    t = lift(mod(g,p)^(c*r))

    M = matrix(3)
    M[0,0] = p
    M[1,0] = inverse_mod(t,p)
    M[1,1] = 1
    M[2,0] = a1 - inverse_mod(t,p)*a2
    M[2,2] = 2^64

    N = M.LLL()

    print("a1=",hex(a1))
    print("a2=",hex(a2))

    return "solution (k1, k2) is given by", (hex(b1),hex(b2))
    #print(N)
    
    #print(mod(b1-inverse_mod(t,p)*b2+a1-inverse_mod(t,p)*a2,p))

# -------------------  Other code: setting parameters and other --------------------

def gen_curve():
    p = 0xffffffffffffffc5
    done = False
    i = 1
    while not done:
        print(i)
        F = FiniteField(p)
        C = EllipticCurve([F(0),F(3)])
        if is_prime(C.cardinality()):
            done = True
            return p
        else:
            p = previous_prime(p)
            i += 1
                
def ecdsa_params():
    p = 0xffffffffffffd21f
    F = FiniteField(p)
    C = EllipticCurve([F(0),F(3)])
    n = 0xfffffffefa23f437
    G = C.lift_x(1)# (1,2)
    x = 0x34aad140ec2c3a3
    return p,F,C,n,G,x

def dsa_params():
    g = 0x17dfdbf2bbbae7d6c052c2fdc5d3288d
    p = 0x89524bfca958c9165a087cc4f889a08f
    q = 0xffffffffffffffc5
    y = 0x2410f15634222d3300eabeb44226cea8
    x = 0x38dbefc062cd4cf3

def dh_params():
    p = 0xffffffffffffffffffffffffffffc3a7
    g = 2
    return p,g

def safe_prime(l=128):
    p = previous_prime(2^l)
    done = False
    i = 0
    while not done:
        print(i)
        if is_prime(Integer((p-1)/2)):
            done = True
            return p
        else:
            p = previous_prime(p)
            i += 1


    
def btoi(b):
    return int.from_bytes(b, "big")

def itob(i, baselen):
    return int.to_bytes(int(i), length=baselen, byteorder="big")

def sign(h, klen=32, return_k=False):
    p,F,C,n,G,x = ecdsa_params()
    d = x
    hi = btoi(h)
    k = ZZ.random_element(2 ** klen)
    r = Integer((G * k).xy()[0])
    s = lift(inverse_mod(k, n) * mod(hi + d * r, n))
    sig = bytes.hex(itob(r, 8)) +" "+ bytes.hex(itob(s, 8))
    if return_k:
        return k,sig
    else:
        return sig

def gen_dsa_prime():
    p = 2*q*random_prime(2^64)+1
    i = 1
    while not is_prime(p):
        p = 2*q*random_prime(2^64)+1
        i += 1
        print(i)
    return p

def gen_sig():
    h = itob(ZZ.random_element(2^64),64/8)
    return bytes.hex(h), sign(h)

#----------------------------Generating plots in the survey ---------------------------------


def kangaroo(lsbs=97):
    p = 65267
    q = 32633
    g = 3
    
    m = 0x1400
    a = 5120+lsbs
    w = 256
    A = lift(mod(g,p)^a)

    c = 4
    S = [(lift(mod(g,p)^s),s) for s in [1,3,7,10]]
    
    a0 = m+w/2
    x0 = lift(mod(g,p)^(a0))
    
    y0 = A
    b0 = 0

    tame_dict = {}

    def H(x):
        return lift(mod(x,c))
    
    xi = x0; ai=a0
    yi = y0; bi=b0

    xmin = 1.5; xmax = 17
    emax = m+w-a-100
    def ct(e):
        return xmin+float(e/(emax)*(xmax-xmin))

    print("\\draw[|-] (%.2f,4.0) -- (%.2f,4);"%(0,1))
    print("\\draw[dotted] (%.2f,4.0) -- (%.2f,4);"%(1,1.5))
    print("\\draw[--] (%.2f,4.0) -- (%.2f,4);"%(1.5,xmax-2))
    print("\\draw[dotted] (%.2f,4.0) -- (%.2f,4);"%(xmax-2,xmax-1))
    print("\\draw[-|] (%.2f,4.0) -- (%.2f,4);"%(xmax-1,xmax))
    print("\\node[anchor=south west] at (%.2f,5) {$m=\mathtt{%s}$};"%(0,hex(m)))
    print("\\node[] at (%.2f,3.5) {$m+w=\mathtt{%s}$};"%(xmax,hex(m+w)))
			
    print("\\node[text=blue] at (%.2f,3.0) {$a$};"%(ct(b0)));
    print("\\node[text=red] at (%.2f,5) {$\mathtt{%s}$};"%(ct(a0-a),str(hex(a0))));
    
    for i in range(4):
        xi = lift(mod(x0*S[H(x0)][0],p))
        ai = lift(mod(a0+S[H(x0)][1],q))
        tame_dict[xi] = ai
        print("\\draw[red,->,bend left] (%.2f,4) to[out=90,in=90] (%.2f,4);"%(ct(a0-a),ct(ai-a)))
        print("\\node[text=red] at (%.2f,5) {$\mathtt{%s}$};"%(ct(ai-a),str(hex(ai))));
        x0 = xi
        a0 = ai
    print("\\draw[red,->,bend left] (%.2f,4) to[out=90,in=90] (%.2f,4);"%(ct(ai-a),ct(ai-a+1)))
    print("\\draw[red,->,bend left] (%.2f,4) to[out=90,in=90] (%.2f,4);"%(ct(ai-a+1),ct(ai-a+2)))
    print("\\draw[red,->,bend left] (%.2f,4) to[out=90,in=90] (%.2f,4);"%(ct(ai-a+2),ct(ai-a+4)))

    for i in range(16):
        yi = lift(mod(y0*S[H(y0)][0],p))
        bi = lift(mod(b0+S[H(y0)][1],q))
        print("\\draw[blue,->,bend left] (%.2f,4) to[out=270,in=270] (%.2f,4);"%(ct(b0),ct(bi)))
        #if i in [1,3,4,5,7,8]:
        print("\\node[text=blue] at (%.2f,3.0) {a+$\mathtt{%s}$};"%(ct(bi),str(hex(bi))));
        if yi in tame_dict:
            ai = tame_dict[yi]
            candidate_a = lift(mod(ai-bi,q))
            print("Success:",i,candidate_a==a,candidate_a,a)
            break


        y0 = yi
        b0 = bi
        
def gen_tree_pq(N=899,p='?11?1',q='?1?0?'):
    candidate_list = []
    candidate_list.push('','',0)
    while candidate_list:
        pass

def gen_children_pq(pi,qi,N=899,p='?11?1',q='?1?0?'):
    i = len(pi)+1
    candidate_list = [(bp+pi,bq+qi) for (bp,bq) in [('0','0'),('0','1'),('1','0'),('1','1')]]
    for new_pi,new_qi in candidate_list:
        if len(new_pi) <= len(p) and p[-i] != '?' and p[-i] != new_pi[-i]:
            continue
        if len(new_qi) <= len(q) and q[-i] != '?' and q[-i] != new_qi[-i]:
            continue
        if mod(N,2^i) != mod(Integer(new_pi,2)*Integer(new_qi,2),2^i):
            new_qi += 'x'
        print('p='+new_pi,'q='+new_qi)
        
def gen_children_dpdq(dpi,dqi,pi,qi,N=899,dp = '?0??1', dq='???0?',kp = 13, kq = 3,e=17):
    i = len(dpi)+1
    print("i=",i)
    dpdq_candidates = [(bp+dpi,bq+dqi) for (bp,bq) in [('0','0'),('0','1'),('1','0'),('1','1')]]
    pq_candidates = [(bp+pi,bq+qi) for (bp,bq) in [('0','0'),('0','1'),('1','0'),('1','1')]]
    for new_dpi,new_dqi in dpdq_candidates:
        if len(new_dpi) <= len(dp) and dp[-i] != '?' and dp[-i] != new_dpi[-i]:
            continue
        if len(new_dqi) <= len(dq) and dq[-i] != '?' and dq[-i] != new_dqi[-i]:
            continue
        for new_pi,new_qi in pq_candidates:
            if mod(e*Integer(new_dpi,2)-1+kp,2^i) != mod(kp*Integer(new_pi,2),2^i):
                continue
            if mod(e*Integer(new_dqi,2)-1+kq,2^i) != mod(kq*Integer(new_qi,2),2^i):
                continue
            if mod(N,2^i) != mod(Integer(new_pi,2)*Integer(new_qi,2),2^i):
                new_qi += 'x'
            print('dp='+new_dpi,'dq='+new_dqi,'p='+new_pi,'q='+new_qi)