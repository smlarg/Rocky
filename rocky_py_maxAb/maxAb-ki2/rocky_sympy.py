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


def floatify(ex):
  if ex.args == ():
    return ex*1.0
  new_args = []
  for arg in ex.args:
    arg = floatify(arg)
    new_args.append(arg)
  return ex.func(*new_args)

ki = 2
kay = sqrt(ki**2 + 1) + ki*I

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

fixRe = lambda func : simplify(de_hyper(re(func.rewrite(sin))))

Ev = fixRe(diff(-Av,t))
Ep = fixRe(diff(-Ap,t))

Bv = fixRe(diff(Av,x))
Bp = fixRe(diff(Ap,x))

print \
'  def vector_potential(x, t):'  + '\n' + \
'    if x < 0:' + '\n' + \
'      return ' + str(simplify(floatify(fixRe(Av)))) + '\n' + \
'    else:' + '\n' + \
'      return ' + str(simplify(floatify(fixRe(Ap))))

print \
'  if (x <= 0) then\n' + \
'    E(2) = ' + str(simplify(floatify(Ev))) + '\n' + \
'  else\n' + \
'    E(2) = ' + str(simplify(floatify(Ep))) + '\n' + \
'  endif'

print \
'  if (x <= 0) then\n' + \
'    B(3) = ' + str(simplify(floatify(Bv))) + '\n' + \
'  else\n' + \
'    B(3) = ' + str(simplify(floatify(Bp))) + '\n' + \
'  endif'


