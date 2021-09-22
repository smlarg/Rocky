This code is an implementation of the Boris particle pusher in Fortran for speed with a python UI for, you know, pythoness.
The fields are not updated by the particles, they're calculated analytically from a given plasma density and absorption.
(Don't try to figure out the latter from collisional plasma parameters, it's...not reasonable; it's meant as a simplified analytic mock-up of hot-particle absorption.)

The field geometry is calculated in python, using sympy, which then writes Fortran code to be compiled; so it's kind of a poor-man's JIT compilation.

Everything remains pretty manual, because you know, graduate student; options for automating it, in either shell or python, are so obvious I won't insult you by explaining or doing them.
