import sys
sys.path.append('.')
import quad

r = quad.Roots()
quad.quadratic(2, 5, 2, r)
print(r.x1, r.x2)
