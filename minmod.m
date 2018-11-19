function c = minmod(a, b)

c = 0.5*(sign(a)+sign(b)).*min(abs(a), abs(b));