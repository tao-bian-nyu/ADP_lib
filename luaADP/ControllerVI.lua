-- Lua part binding, load necessary interface.
ffi = require('ffi')
C = ffi.C
ffi.cdef[[
typedef struct ControllerVI ControllerVI;
typedef struct SymmetricMatrix SymmetricMatrix;
typedef struct Step Step;

ControllerVI* ControllerVI_new();

double* ControllerVI_learner(ControllerVI* self, const double* x, const double* u, const int n, const int m, const double dt, const double t);

void ControllerVI_dispAll(ControllerVI* self);

void ControllerVI_delete(ControllerVI* self);
void free(double*);
]]
--controllerVI = ffi.load('libControllerVI.so')
controllerVI = ffi.load('luaADP/libControllerVI.so')

controllerVI_index = {
	--Input = controllerVI.ControllerVI_input,
	--Learner = controllerVI.ControllerVI_learner
}
--controllerVI_mt = ffi.metatype('ControllerVI', {
	--__index = controllerVI_index
--})

controllerVI_index.__index = controllerVI_index

function ControllerVI(...)
	--ControllerVI = ffi.gc(controllerVI.ControllerVI_new(), controllerVI.ControllerVI_delete)
	local self = {super = controllerVI.ControllerVI_new(...)}
	ffi.gc(self.super, controllerVI.ControllerVI_delete)
	return setmetatable(self, controllerVI_index)
	--return ControllerVI
end

function controllerVI_index.learn(self,x,u,n,m,dt,t)
	--local name = controllerVI.ControllerVI_learner(...)
	super2 = controllerVI.ControllerVI_learner(self.super,x,u,n,m,dt,t)
	ffi.gc(super2, C.free)
	return super2
end
