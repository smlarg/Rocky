from sympy import *

def de_hyper(ex):
  if ex.args == ():
    return ex
  new_args = []
  for arg in ex.args:
    arg = de_hyper(arg)
    if arg.func in (sinh, cosh):
      new_args.append(arg.rewrite(exp, deep = False))
    else:
      new_args.append(arg)
  return ex.func(*new_args)


def realify(ex):
  if ex.args == ():
    return ex*1.0
  new_args = []
  for arg in ex.args:
    arg = realify(arg)
    new_args.append(arg)
  return ex.func(*new_args)


kay = Rational(1,2) + I

k = symbols('k')
x, t = symbols('x t', real = True)
Ai = 6
ar = (1-k)/(1+k)
at = 2/(1+k)
ar = ar.subs(k, kay)
at = at.subs(k, kay)
Av = Ai * (exp(I*(x-t)) + ar * exp(I*(-x-t)))
Ap = Ai * at * exp(I*(k*x - t))

Ap = Ap.subs(k, kay)
Ap = Ap.rewrite(sin)
Ap = re(Ap)
Ap = de_hyper(Ap)

Av = simplify(re(Av))
Ap = simplify(re(Ap))

Ev = simplify(re(diff(-Av,t)))
Ep = simplify(re(diff(-Ap,t)))

Bv = simplify(re(diff(Av, x)))
Bp = simplify(re(diff(Ap, x)))

print \
'  def vector_potential(x, t):' + '\n' + \
'    if x < 0:' + '\n' + \
'      return ' + str(simplify(realify(Av))) + '\n' + \
'    else:' + '\n' + \
'      return ' + str(simplify(realify(Ap)))

print \
'  if (x <= 0) then' + '\n' + \
'    E(2) = ' + str(simplify(realify(Ev))) + '\n' + \
'  else' + '\n' + \
'    E(2) = ' + str(simplify(realify(Ep))) + '\n' + \
'  endif'

print \
'  if (x <= 0) then' + '\n' + \
'    B(3) = ' + str(simplify(realify(Bv))) + '\n' + \
'  else' + '\n' + \
'    B(3) = ' + str(simplify(realify(Bp))) + '\n' + \
'  endif'


