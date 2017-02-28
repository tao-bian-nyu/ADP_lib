controller= require('ControllerVI')
ffi = require('ffi')
matrix = require('matrix')
control= ControllerVI()
c = control:learn(ffi.new("double[3]", {1,2,3}),ffi.new("double[1]", {1}), 3, 1, 1,2)


A = matrix{{-1,2,1},{0.1,-1,2},{0,0.1,1}}
B = matrix{{0},{0},{1}}
x = matrix{{2},{0},{0}}
K = matrix{{2,0,0}}

matrix.print(B*2)
dt = 0.0001
t = 0

for i = 1,10 do
	u = -matrix.mul(K,x)
	dx = (matrix.mul(A,x) + matrix.mul(B,u)) * dt
	x = x + dx
	K = control:learn(ffi.new("double[3]", x),ffi.new("double[1]", u), 3, 1, dt,t)
	i = i+1
	t = t + dt
end
--c = learn(control,ffi.new("double[3]", {1,2,3}),ffi.new("double[1]", {1}), 3, 1, 1,2)
--print(c[1])

