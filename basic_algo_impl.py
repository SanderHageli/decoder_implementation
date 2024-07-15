from sage.all import *
import random

q = Integer(4)
r = Integer(37)
d = r-10
g = Integer(int((q**2-q)/2))
k.<z> = GF(q**2,"z", modulus="primitive", repr="log")
F = VectorSpace(k,q**3)

def gen_riemann_roch_basis(m):
    """Generates a basis for the Riemann-Roch space L(D,mQ)"""
    basis, i, j = [], 0, 0
    while q*(i-j) + (q+1)*j <= m:
        for j in range(i+1):
            if j<q and q*(i-j) + (q+1)*j <= m:
                basis.append(lambda x,y, exp1=i-j, exp2=j: (x**exp1)*(y**exp2)) #basis is generated is such a way to get the order 1,x,y,x^2,xy,y^2,...
        j,i=0,i+1
    return basis

def find_curve_points():
    points = []
    for i in k:
        for j in k:
            if j**q + j == i**(q+1): #cycles through the field elements to see if they satisfies the equation
                points.append((i,j))
    return points

def generate_error(m):
    """Generates a vector with random field elements at m random vector indices"""
    error_indices, error_vector = random.sample(range(q**3),m), []
    for  i in range(q**3):
        if i in error_indices:
            error_vector.append(k.random_element())
        else:
            error_vector.append(0)
    return error_indices,F(error_vector)

def parity_matrix():
    M, v = MatrixSpace(k,len(riemann_roch_basis),q**3), []
    for h in c_r_basis:
        v+=h
    return M(v)

def evaluate(f, points):
    return F([f(point[0],point[1]) for point in points])

def generate_codeword():
    riemann_roch_basis = gen_riemann_roch_basis(q**3+q**2-q-2-r)
    coeffs = [k.random_element() for i in riemann_roch_basis]
    g = lambda x,y: sum([coeffs[i]*riemann_roch_basis[i](x,y) for i in range(len(coeffs))]) #inner product compatible with the evaluation funciton
    return evaluate(g,points)

def syndrome(i,j,word):
    return word.inner_product(c_r_basis[i].pairwise_product(c_r_basis[j]))

def syndrome_matrix(a,b,word):
    M, v = MatrixSpace(k,a,b), []
    for i in range(a):
        v+=[syndrome(i,j,word) for j in range(b)]
    return M(v)

def N(f):
    """Finds the correpsonding indecies where f is zero at"""
    g = lambda x,y: sum([f[i]*riemann_roch_basis[i](x,y) for i in range(len(f))])
    g_evaluate, zero_indices = evaluate(g,points), []
    for i in range(q**3):
        if g_evaluate[i]==0:
            zero_indices.append(i)
    return zero_indices

def solve_matrix_system(y,H,indices: list):
    """Solves the linear system xH=yH with y being zero at the indices spesified in indices"""
    HT,v = H.transpose(), []
    for i in indices:
        v+= [j for j in HT[i]]
    reduced_matrix = MatrixSpace(k,len(indices),len(HT[0]))(v).transpose()
    solution = [i for i in reduced_matrix.solve_right(H*y)]
    
    append_zeroes = [0]*(q**3)
    for i in range(q**3):
        if i in indices:
            append_zeroes[i] = solution[0]
            del solution[0]
    return F(append_zeroes)

def decoder(word):
    ker = syndrome_matrix(r_1-g+1, r-r_1-g+1, word).right_kernel()
    if len(ker)==1:
        return "Cannot decode. Kernel is trivial"
    elif len(N(ker[1]))>=d:
        return "Cannot decode. N(f) has more than d elements"
    else:
        H = parity_matrix()
        potential_error = solve_matrix_system(word,H,N(ker[1]))
        return word-potential_error

points = find_curve_points()
riemann_roch_basis = gen_riemann_roch_basis(r)
c_r_basis = [evaluate(f, points) for f in riemann_roch_basis]

n = Integer(int((d-1-g)/2))
r_1 = 21

codeword = generate_codeword()
error_indices,error_vector = generate_error(n)
y = codeword + error_vector

x = decoder(y)
print(x==codeword)