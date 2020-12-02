import numpy as np
from scipy.linalg import solve, norm, det, eigvals,qr, eig
# from scipy.optimize import newton

# 予測子
def predict(y, h, fx, fa, w):
    n = len(w)
    A = np.zeros((n, n))
    A[:-1, :-1] = fx(y[:-1], np.array([y[-1]]))
    A[:-1, -1] = fa(y[:-1], np.array([y[-1]]))
    A[-1,:] = w
    b = np.zeros(n)
    b[-1] = 1.0
    v = solve(A, b)
    v /= norm(v)
    return y+h*v, v


# 修正子
def correct(x, h, y, v, f, fx, fa, eps=1.0e-8, maxiter=10):
    def F(x, h, y, v, f, fx, fa):
        ans = np.zeros(len(x))
        ans[:-1] = f(x[:-1],np.array([x[-1]]))
        ans[-1] = np.dot(x - y, v) - h
        return ans

    def DF(x, h, y, v, f, fx, fa):
        A = np.zeros((len(x), len(x)))
        A[:-1, :-1] = fx(x[:-1],np.array([x[-1]]))
        A[:-1, -1] = fa(x[:-1],np.array([x[-1]]))
        A[-1,:] = v
        return A
    
    z = np.copy(x)
    b = F(z, h, y, v, f, fx, fa)
    for m in range(maxiter):
        A = DF(z, h, y, v, f, fx, fa)
        z -= solve(A, b)
        b = F(z, h, y, v, f, fx, fa)
        if norm(b) < eps:
            return z
    print("not converged")


# 解の追跡バージョン0
def pathfollow00(x, a, func, dfdx, dfda, nmax=100, h=0.1, epsr=1.0e-10, epsb=1.0e-10):
    bd = []
    ndim = len(x) + len(a)
    y = np.zeros(ndim)
    y[:-1] = x
    y[-1] = a[0]
    w = np.zeros(ndim)
    w[-1] = 1.0
    x0 = y[:-1]
    a0 = np.array([y[-1]])
    # 予測子
    yt,v=predict(y, h, dfdx, dfda, w)
    bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v})
    print("# TY", "a", "x")

    for m in range(1, nmax+1):
        y = correct(yt, h, y, v, func, dfdx, dfda, eps=epsr)
        w = v
        x0 = y[:-1]
        a0 = np.array([y[-1]])
        # 予測子
        yt, v = predict(y, h, dfdx, dfda, w)

        bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v})
        print("R", a0[0], x0[0])
    return bd


# 極限点の特定
def locateLP(y, v, f, fx, fa, eps=1.0e-8, h=1.0e-8, maxiter=10):
    n = len(y)
    z = np.copy(y)
    e = np.zeros(n)
    e[-1] = 1.0
    d = fa(z[:-1], np.array([z[-1]]))
    c = v[:-1]
    for m in range(maxiter):
        x0 = z[:-1]
        a0 = np.array([z[-1]])
        b = np.zeros(len(y))
        b[:-1] = f(x0, a0)
        M = np.zeros((n, n))
        M[:-1,:-1] = fx(x0, a0)
        M[-1,:-1] = c
        M[:-1,-1] = d
        w = solve(M, e)
        b[-1] = w[-1]
        if norm(b) < eps:
            return z
        phi = solve(M.T, e)
        A = np.zeros((n, n))
        A[:-1,:-1] = M[:-1, :-1]
        A[:-1,-1] = fa(x0, a0)
        for j in range(n-1):
            ej = np.zeros(n-1)
            ej[j] = 1.0
            M[-1, j] = -phi[:-1] @ ((fx(x0+h*ej, a0) - M[:-1,:-1]) / h) @ w[:-1]
        M[-1,-1] = -phi[:-1] @ ((fx(x0, a0+h) -  M[:-1,:-1]) / h) @ w[:-1]
        z -= solve(M, b)
    print("not converged")


# 分岐点の特定
def locateBP(y, v, f, fx, fa, eps=1.0e-8, maxiter=10, h=1.0e-7):
    n = len(y)-1
    z = np.zeros(2*(n+1))
    z[:n] = y[:-1]
    z[n] = y[-1]
    z[n+1] = 0.0
    # QR factorization with pivotting
    Q, R, P = qr(fx(z[:n], np.array([z[n]])), mode='full', pivoting=True)
    en = np.zeros(n)
    en[-1] = 1.0
    z[n+2:] = Q @ en

    for m in range(maxiter):
        x0 = z[:n]
        a0 = np.array([z[n]])
        beta0 = np.array([z[n+1]])
        p0 = z[n+2:]
        b = np.zeros(len(z))
        b[:n] = f(x0, a0) + beta0 * p0
        A = fx(x0, a0)
        c = fa(x0, a0)
        b[n:2*n] = A.T @ p0
        b[2*n] = p0 @ c.T
        b[2*n+1] = p0 @ p0 - 1.0
        if norm(b) < eps:
            break
        M = np.zeros((2*(n+1), 2*(n+1)))
        M[:n,:n], M[:n, n], M[:n,n+1]  = A, c, p0
        for j in range(n):
            M[j, n+2+j] = beta0
            ej = np.zeros(n)
            ej[j] = 1.0
            M[n:2*n, j] = np.transpose((fx(x0+h*ej, a0) - A) / h) @ p0
            M[2*n, j] = p0 @ np.transpose(fa(x0+h*ej, a0) - c) / h
        M[n:2*n, n] = np.transpose((fx(x0, a0+h) - A) / h) @ p0
        M[n:2*n, n+2:] = A.T
        M[2*n, n] = p0 @ np.transpose(fa(x0, a0+h) - c) / h
        M[2*n, n+2:] = c
        M[2*n+1, n+2:] = 2 * p0
        z -= solve(M, b)
    if m == maxiter:
        print("not converged")
    else:
        phi = np.zeros(n+1)
        phi[:-1] = z[n+2:]
        return z[:n+1], phi


# 分岐点における接ベクトル計算
def calcTangentVectorBP(y, v0, v1, phi, f, eps=1.0e-8):
    h = 1.0e-6
    a = 0.0
    b = 1.0
    c = 0.5
    Ba = phi[:-1] @ (f(y[:1] + h*v1[:-1], np.array([y[-1]+h*v1[-1]])) + f(y[:1] - h*v1[:-1], np.array([y[-1]-h*v1[-1]]))) / h**2
    Bb = phi[:-1] @ (f(y[:1] + h*v1[:-1], np.array([y[-1]+h*v1[-1]])) + f(y[:1] - h*v1[:-1], np.array([y[-1]-h*v1[-1]]))) / h**2
    while abs(b-a) < eps:
        c = 0.5 * (a + b)
        xc = y[:-1]
        ac=np.array([y[-1]])
        vcx = (1-c)*v0[:-1] + c * v1[:-1]
        vca = np.array([(1-c)*v0[-1] + c * v1[-1]])
        Bc = phi[:-1] @ (f(xc+h*vcx, ac+h*vca) + f(xc-h*vcx, ac-h*vca)) / h**2
        if Ba * Bc > 0:
            a = c
            Ba = Bc
        elif Bc * Bb > 0:
            b = c
            Bb = Bc
        else:
            return (1-c)*v0+c*v1
    return (1-c)*v0+c*v1   


# 分岐点における追跡枝切替方向の計算
def calcSwitchingVectorBP(pt, f, fx, fa, eps=1.0e-8, problem=None, period=1):
    if pt['TY'] != 'B':
        print("Not a bifurcation point")

    if problem == 'map':
        def func(x, a):
            y = np.copy(x)
            for m in range(0,period):
                y = f(y,a)
            return y - x
        def dfdx(x,a):
            A = fx(x,a)
            y = np.copy(x)
            for m in range(1,period):
                y = f(y, a)
                A = fx(y, a) @ A
            return A - np.identity(len(x))
        def dfda(x,a):
            b = fa(x,a)
            y = np.copy(x)
            for m in range(period):
                y = f(y, a)
                b = fx(y, a) @ b + fa(y, a)
            return b
    else:
        func = f
        dfdx = fx
        dfda = fa

    h = 1.0e-6
    x = np.copy(pt['x'])
    a = np.copy(pt['a'])
    v = np.copy(pt['v'])
    phi = np.copy(pt['phi'])
    n = len(x)+1
    I = np.identity(n)
    A = np.zeros((n, n))
    A[:-1,:-1] = dfdx(x,a)
    A[:-1,-1] = dfda(x,a)
    A[-1,:] = v
    ## QR
    Q,R,P = qr(A.T, mode='full', pivoting=True)
    en = np.zeros(n)
    en[-1] = 1.0
    z = np.zeros(n+1)
    z[:-1] = Q @ en
    z[-1] = 0.0
    b = np.zeros(n+1)
    b[:-1] = (A - z[-1]*I) @ z[:-1]
    b[-1] = 0.5 * (z[:-1]@z[:-1] - 1.0)
    for m in range(15):
        B = np.zeros((n+1,n+1))
        B[:-1,:-1] = A - z[-1] * I
        B[:-1,-1] = -z[:-1]
        B[-1,:-1] = z[:-1]
        z -= solve(B, b)
        b[:-1] = (A - z[-1]*I) @ z[:-1]
        b[-1] = 0.5 * (z[:-1]@z[:-1] - 1.0)
        if norm(b) < eps:
            break
    if m == 15:
        print("not converged")
    v2 = np.copy(z[:-1])
    a22 = phi[:-1] @ (func(x+h*v[:-1], a+h*v[-1]) + func(x-h*v[:-1], a-h*v[-1])) / h**2
    Bp = (func(x+h*(v[:-1]+v2[:-1]), a+h*(v[-1]+v2[-1])) + func(x-h*(v[:-1]+v2[:-1]), a-h*(v[-1]+v2[-1]))) / h**2
    Bm = (func(x+h*(v[:-1]-v2[:-1]), a+h*(v[-1]-v2[-1])) + func(x-h*(v[:-1]-v2[:-1]), a-h*(v[-1]-v2[-1]))) / h**2
    a12 = phi[:-1] @ (Bp - Bm) / 4
    w = -0.5*a22*v/a12 + v2 
    return w/norm(w)


# 周期倍分岐点の特定
def locatePD(y, f, fx, fa, period, eps=1.0e-8, h=1.0e-8, maxiter=10):
    def F(x, a):
        z = np.copy(x)
        for m in range(period):
            z = f(z,a)
        return z
    def Fx(x,a):
        A = fx(x,a)
        y = np.copy(x)
        for m in range(1,period):
            y = f(y, a)
            A = fx(y, a) @ A
        return A
    def Fa(x,a):
        b = fa(x,a)
        y = np.copy(x)
        for m in range(1,period):
            y = f(y, a)
            b = fx(y, a) @ b + fa(y, a)
        return b

    n = len(y)-1
    z = np.zeros(2*n+1)
    z[:n+1] = np.copy(y)
    if n > 1:
        Q,R,P = qr(Fx(y[:-1], np.array([y[-1]]))+np.identity(n), pivoting=True)
    else:
        Q = np.array([[1.0]])
    en = np.zeros(n)
    en[-1] = 1.0
    z[n+1:] = Q @ en
    for m in range(maxiter):
        x0 = z[:n]
        a0 = np.array([z[n]])
        b = np.zeros(2*n+1)
        b[:n] = F(x0, a0) - x0
        A = Fx(x0, a0)
        b[n:-1] = (A + np.identity(n)) @ z[n+1:]
        b[-1] = 0.5 * (z[n+1:] @ z[n+1:] - 1.0)

        M = np.zeros((2*n+1, 2*n+1))
        M[:n,:n] = A - np.identity(n)
        M[:n, n] = Fa(x0, a0)
        M[n:-1, n+1:] = A + np.identity(n)
        M[-1, n+1:] = z[n+1:]
        for j in range(n):
            ej = np.zeros(n)
            ej[j] = 1.0
            M[n:-1, j] = ((Fx(x0+h*ej, a0) - A) / h) @ z[n+1:]
        M[n:-1, n] = ((Fx(x0, a0+h) - A) / h) @ z[n+1:]
        z -= solve(M, b)
        if norm(b) < eps:
            break
    if m == maxiter:
        print("not converged")
    else:
        x0, a0 = z[:n], np.array([z[n]])
        Q, R, P = qr(Fx(F(x0, a0), a0) @ Fx(x0, a0) - np.identity(n), mode='full', pivoting=True)
        en = np.zeros(n)
        en[-1] = 1.0
        phi = np.zeros(n+1)
        phi[:-1] = Q @ en
        return z[:n+1], phi


# 周期倍分岐点における追跡枝切替方向の計算
def calcSwitchingVectorPD(pt, f, fx, fa, eps=1.0e-8, period=1):
    if pt['TY'] != 'P':
        print("Not a bifurcation point")

    def func(x, a):
        y = np.copy(x)
        for m in range(0,period):
            y = f(y,a)
        return y - x
    def dfdx(x,a):
        A = fx(x,a)
        y = np.copy(x)
        for m in range(1,period):
            y = f(y, a)
            A = fx(y, a) @ A
        return A - np.identity(len(x))
    def dfda(x,a):
        b = fa(x,a)
        y = np.copy(x)
        for m in range(1,period):
            y = f(y, a)
            b = fx(y, a) @ b + fa(y, a)
        return b

    h = 1.0e-6
    x = np.copy(pt['x'])
    a = np.copy(pt['a'])
    v = np.copy(pt['v'])
    phi = np.copy(pt['phi'])
    n = len(x)+1
    I = np.identity(n)
    A = np.zeros((n, n))
    A[:-1,:-1] = dfdx(x,a)
    A[:-1,-1] = dfda(x,a)
    A[-1,:] = v
    Q,R,P = qr(A.T, mode='full', pivoting=True)
    en = np.zeros(n)
    en[-1] = 1.0
    z = np.zeros(n+1)
    z[:-1] = Q @ en
    z[-1] = 0.0
    b = np.zeros(n+1)
    b[:-1] = (A - z[-1]*I) @ z[:-1]
    b[-1] = 0.5 * (z[:-1]@z[:-1] - 1.0)
    for m in range(15):
        B = np.zeros((n+1,n+1))
        B[:-1,:-1] = A - z[-1] * I
        B[:-1,-1] = -z[:-1]
        B[-1,:-1] = z[:-1]
        z -= solve(B, b)
        b[:-1] = (A - z[-1]*I) @ z[:-1]
        b[-1] = 0.5 * (z[:-1]@z[:-1] - 1.0)
        if norm(b) < eps:
            break
    if m == 15:
        print("not converged")
    v2 = np.copy(z[:-1])
    a22 = phi[:-1] @ (func(x+h*v[:-1], a+h*v[-1]) + func(x-h*v[:-1], a-h*v[-1])) / h**2
    Bp = (func(x+h*(v[:-1]+v2[:-1]), a+h*(v[-1]+v2[-1])) + func(x-h*(v[:-1]+v2[:-1]), a-h*(v[-1]+v2[-1]))) / h**2
    Bm = (func(x+h*(v[:-1]-v2[:-1]), a+h*(v[-1]-v2[-1])) + func(x-h*(v[:-1]-v2[:-1]), a-h*(v[-1]-v2[-1]))) / h**2
    a12 = phi[:-1] @ (Bp - Bm) / 4
    w = -0.5*a22*v/a12 + v2 
    return w/norm(w)


# ホップ分岐点の特定
def locateHB(y, v, w, omega, func, dfdx, dfda, eps, h=1.0e-8):
    n = len(y) - 1
    z = np.zeros(3*n+2)
    I = np.identity(n)
    z[:n+1] = np.copy(y)
    z[n+1:2*n+1] = np.copy(v)
    z[2*n+1:3*n+1] = np.copy(w)
    z[-1] = omega
    for m in range(10):
        x0 = z[:n]
        a0 = np.array([z[n]])
        v0 = z[n+1:2*n+1]
        w0 = z[2*n+1:3*n+1]
        omg0 = z[-1]
        f = np.zeros(3*n+2)
        f[:n] = func(x0, a0)
        A = dfdx(x0, a0)
        f[n:2*n] = A @ v0 + omg0 * w0
        f[2*n:3*n] = A @ w0 - omg0 * v0
        f[-2] = v0 @ v0 + w0 @ w0 - 1
        f[-1] = v0 @ w0
        J = np.zeros((3*n+2, 3*n+2))
        J[:n, :n] = A
        J[:n, n] = dfda(x0, a0)

        J[n:2*n, n+1:2*n+1] = A
        J[n:2*n, 2*n+1:3*n+1] = omg0 * I
        J[n:2*n, -1] = w0

        for j in range(n):
            ej = np.zeros(n)
            ej[j] = 1.0
            J[n:2*n, j] = ((dfdx(x0+h*ej, a0) - A) / h) @ v0
            J[2*n:3*n, j] = ((dfdx(x0+h*ej, a0) - A) / h) @ w0
        J[n:2*n, n] = ((dfdx(x0, a0+h) - A) / h) @ v0
        J[2*n:3*n, n] = ((dfdx(x0, a0+h) - A) / h) @ w0

        J[2*n:3*n, n+1:2*n+1] = -omg0 * I
        J[2*n:3*n, 2*n+1:3*n+1] = A
        J[2*n:3*n, -1] = -v0

        J[-2, n+1:2*n+1] = 2 * v0
        J[-2, 2*n+1:3*n+1] = 2 * w0
        J[-1, n+1:2*n+1] = w0
        J[-1, 2*n+1:3*n+1] = v0

        z -= solve(J, f)
        if norm(f) < eps:
            return z[:n+1], z[n+1:2*n+1], z[2*n+1:3*n+1],z[-1] 
    print("locateHB: not converged")



# 解の追跡（以上をすべて盛り込んだもの）
def pathfollow(x, a, f, fx, fa, w=None,nmax=100, h=0.1, epsr=1.0e-10, epsb=1.0e-10, amin=-np.inf, amax=np.inf, problem=None, period=1, quiet=None):
    if problem == 'map':
        def func(x, a):
            z = np.copy(x)
            for m in range(period):
                z = f(z,a)
            return z - x
        def dfdx(x,a):
            A = fx(x,a)
            y = np.copy(x)
            for m in range(1,period):
                y = f(y, a)
                A = fx(y, a) @ A
            return A - np.identity(len(x))
        def dfda(x,a):
            b = fa(x,a)
            y = np.copy(x)
            for m in range(1,period):
                y = f(y, a)
                b = fx(y, a) @ b + fa(y, a)
            return b
    else:
        func = f
        dfdx = fx
        dfda = fa

    bd = []
    ndim = len(x) + len(a)
    y = np.zeros(ndim)
    y[:-1] = x
    y[-1] = a[0]
    x0 = y[:-1]
    a0 = np.array([y[-1]])
    # 極限点のテスト関数
    flp0 = det(dfdx(x0, a0))
    # 最初の予測
    if w is None:
        w = np.zeros(ndim)
        w[-1] = 1.0
        yt,v=predict(y, h, dfdx, dfda, w)
    else:
        v = np.copy(w)
        yt = y + h * w
    # 分岐点のテスト関数
    B = np.zeros((ndim, ndim))
    B[:-1,:-1] = dfdx(x0,a0)
    B[:-1,-1]=dfda(x0,a0)
    B[-1,:]=v
    fbp0 = det(B)
    # 周期倍分岐のテスト関数
    if problem == 'map':
        fpd0 = det(dfdx(x0,a0)+2*np.identity(len(x)))
    # 安定性
    if problem == 'equilibria':
        ev, vr = eig(dfdx(x0, a0))
        udim = np.count_nonzero(np.real(ev) > 0.0)
        bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v, 'udim':udim})
        # Hopf分岐のテスト関数
        fhb0 = np.product([np.product([ev[i]+ev[j] for i in range(j+1,len(ev))]) for j in range(len(ev)-1)])
    elif problem == 'map':
        ev = norm(eigvals(fx(x0, a0)), ord=np.inf)
        udim = np.count_nonzero(ev > 1.0)
        bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v, 'udim':udim})
    else:
        bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v})

    if not quiet:
        print("# TY", "a", "x")

    for m in range(1, nmax+1):
        y = correct(yt, h, y, v, func, dfdx, dfda, eps=epsr)
        w = v
        x0 = y[:-1]
        a0 = np.array([y[-1]])
        # 極限点のテスト関数
        flp1 = det(dfdx(x0, a0))
        # 予測子
        yt, v = predict(y, h, dfdx, dfda, w)
        # 分岐点のテスト関数
        B = np.zeros((ndim, ndim))
        B[:-1,:-1] = dfdx(x0,a0)
        B[:-1,-1]=dfda(x0,a0)
        B[-1,:]=v
        fbp1 = det(B)
        # 周期倍分岐点のテスト関数
        if problem == 'equilibria':
            ev, vr = eig(dfdx(x0, a0))
            # Hopf分岐のテスト関数
            fhb1 = np.product([np.product([ev[i]+ev[j] for i in range(j+1,len(ev))]) for j in range(len(ev)-1)])
        elif problem == 'map':
            fpd1 = det(dfdx(x0,a0)+2*np.identity(len(x)))

        if m > 1 and fbp0*fbp1 < 0:
            ybp,phi = locateBP(y, w, func, dfdx, dfda, eps=epsb)
            vbp = calcTangentVectorBP(ybp, w, v, phi, func, eps=1.0e-10)
            if problem == 'equilibria':
                ev_bp = eigvals(dfdx(ybp[:-1], np.array([ybp[-1]])))
                udim = np.count_nonzero(np.real(ev_bp) > 0.0)
                bd.append({'TY':"B", 'x':ybp[:-1], 'a':np.array([ybp[-1]]), 'v':vbp, 'phi':phi, 'udim':udim})
            elif problem == 'map':
                ev_bp = norm(eigvals(fx(ybp[:-1], np.array([ybp[-1]]))), ord=np.inf)
                udim = np.count_nonzero(ev_bp > 1.0)
                bd.append({'TY':"B", 'x':ybp[:-1], 'a':np.array([ybp[-1]]), 'v':vbp, 'phi':phi, 'udim':udim})
            else:
                bd.append({'TY':"B", 'x':ybp[:-1], 'a':np.array([ybp[-1]]), 'v':vbp, 'phi':phi})

            if not quiet:
                print("B", ybp[-1], ybp[0])

        elif m > 1 and flp0 * flp1 < 0:
            ylp = locateLP(y, w, func, dfdx, dfda, eps=epsb)
            if problem == 'equilibria':
                ev_lp = np.real(eigvals(dfdx(ylp[:-1], np.array([ylp[-1]]))))
                udim = np.count_nonzero(ev_lp > 0.0)
                bd.append({'TY':"L", 'x':ylp[:-1], 'a':np.array([ylp[-1]]), 'v':v, 'udim': udim})
            elif problem == 'map':
                ev_lp = norm(eigvals(dfdx(ylp[:-1], np.array([ylp[-1]]))), ord=np.inf)
                udim = np.count_nonzero(ev_lp > 1.0)
                bd.append({'TY':"L", 'x':ylp[:-1], 'a':np.array([ylp[-1]]), 'v':v, 'udim': udim})
            else:
                bd.append({'TY':"L", 'x':ylp[:-1], 'a':np.array([ylp[-1]]), 'v':v})

            if not quiet:
                print("L", ylp[-1], ylp[0])

        elif problem == 'equilibria' and m > 1 and fhb0 * fhb1 < 0:
            re_ev = np.abs(np.real(ev))
            i_hb = np.argmin(re_ev)
            ev_hb = ev[i_hb]
            vr_hb = vr[i_hb]
            if np.abs(np.imag(ev_hb)) > 1.0e-15:
                v_hb = np.real(vr[:, i_hb])
                i_v_hb = np.imag(vr[:, i_hb]) if np.imag(ev_hb) > 0 else -np.imag(vr[:, i_hb])
                w_hb = i_v_hb - (v_hb @ i_v_hb) * v_hb / (v_hb @ v_hb)
                vw_scale = np.sqrt(v_hb @ v_hb + w_hb @ w_hb)
                v_hb /= vw_scale
                w_hb /= vw_scale
                omega0 = np.abs(np.imag(ev_hb))
                yhb, vhb, whb, omg = locateHB(y, v_hb, w_hb, omega0, func, dfdx, dfda, eps=epsb)
                bd.append({'TY':"H", 'x':yhb[:-1], 'a':np.array([yhb[-1]]), 'v':v, 'udim': udim, 'vr':vhb, 'vi':whb, 'omg':omg})
                if not quiet:
                    print("H", yhb[-1], yhb[0])
            else:
                pass

        elif problem == 'map' and m > 1 and fpd0 * fpd1 < 0:
            ypd, phi = locatePD(y, f, fx, fa, period, eps=epsb)
            ev_pd = norm(eigvals(dfdx(ypd[:-1], np.array([ypd[-1]]))), ord=np.inf)
            udim = np.count_nonzero(ev_pd > 1.0)
            bd.append({'TY':"P", 'x':ypd[:-1], 'a':np.array([ypd[-1]]), 'v':v, 'udim': udim, 'phi':phi})
            if not quiet:
                print("P", ypd[-1], ypd[0])


        if problem == 'equilibria':
            udim = np.count_nonzero(ev > 0.0)
            bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v, 'udim':udim})
        elif problem == 'map':
            ev = norm(eigvals(dfdx(x0, a0)), ord=np.inf)
            udim = np.count_nonzero(ev > 1.0)
            bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v, 'udim':udim})
        else:
            bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v})

        flp0 = flp1
        fbp0 = fbp1
        if problem == 'equilibria':
            fhb0 = fhb1
        elif problem == 'map':
            fpd0 = fpd1

        if not quiet:
            print("R", a0[0], x0[0], flp1, fbp1)

        if a0 < amin or a0 > amax:
            print("# parameter arrived at boundary")
            break

    # 分岐点と極限点のインデックスをリストに格納する
    if problem == 'map':
        bp = []
        lp = []
        pd = []
        for i,pt in enumerate(bd):
            if (pt['TY'] == 'B'):
                bp.append(i)
            elif (pt['TY'] == 'L'):
                lp.append(i)
            elif (pt['TY'] == 'P'):
                pd.append(i)
        return bd,bp,lp,pd
    elif problem == 'equilibria':
        bp = []
        lp = []
        hb = []
        for i,pt in enumerate(bd):
            if (pt['TY'] == 'B'):
                bp.append(i)
            elif (pt['TY'] == 'L'):
                lp.append(i)
            elif (pt['TY'] == 'H'):
                hb.append(i)
        return bd,bp,lp,hb
    else:
        bp = []
        lp = []
        for i,pt in enumerate(bd):
            if (pt['TY'] == 'B'):
                bp.append(i)
            elif (pt['TY'] == 'L'):
                lp.append(i)
        return bd,bp,lp
