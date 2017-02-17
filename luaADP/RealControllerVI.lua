
-- Lua part binding, load necessary interface.

ffi = require('ffi')
ffi.cdef[[
typedef struct ControllerVI ControllerVI;
typedef struct SymmetricMatrix SymmetricMatrix;
typedef struct Step Step;

ControllerVI* ControllerVI_new(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, Step* stepf);

const double* ControllerVI_learner(ControllerVI* self, const double* x, const double* u, const double dt, const double t);

void ControllerVI_dispAll(ControllerVI* self);

void ControllerVI_delete(ControllerVI* self);
]]
controllerVI = ffi.load('ControllerVI.so')

controllerVI_index = {
	--Input = controllerVI.ControllerVI_input,
	Learner = controllerVI.ControllerVI_learner
}
controllerVI_mt = ffi.metatype('ControllerVI', {
	__index = controllerVI_index
})
function ControllerVI(...)
	--ControllerVI = ffi.gc(controllerVI.ControllerVI_new(), controllerVI.ControllerVI_delete)
	local self = {super = controllerVI.ControllerVI_new(...)}
	ControllerVI = ffi.gc(self.super(...), controllerVI.ControllerVI_delete)
	return ControllerVI
end

