from sympy.geometry import Point3D
from sympy import *
from sympy import symbols, sqrt
from sympy import symbols, Matrix, sqrt, acos, asin, init_printing
from sympy.vector import CoordSys3D, Vector as SympyVector
from IPython.display import Image, display,Math
# Símbolos unitarios
iv = symbols(r'\hat{\mathbf{i}}')
jv = symbols(r'\hat{\mathbf{j}}')
kv = symbols(r'\hat{\mathbf{k}}')
stri='î'
strj='ĵ'
strk='k̂'
from lib_Mathematica import *
from lib_Mathbasic  import *
class MyVector:
    """
    A comprehensive 3D vector class for physics, engineering, and mathematics.
    
    FEATURES:
    ---------
    • Multiple construction methods (components, points, forces, symbolic)
    • Full vector algebra (addition, subtraction, scaling, dot/cross products)
    • Geometric operations (projections, angles, rotations, direction cosines)
    • Physics operations (torque, moments about axes, force calculations)
    • Symbolic computation support via SymPy
    • Unit vector management and coordinate system integration
    
    CONSTRUCTION OPTIONS:
    --------------------
    1. From components: MyVector(1, 2, 3) or MyVector([1, 2, 3])
    2. From two points: MyVector(Point3D(0,0,0), Point3D(1,2,3))
    3. As a force: MyVector([1,1,0], F=10, origin=(1,2,3))
    4. With direction: MyVector(dire=(1,0,0), F=5)
    5. From SymPy objects: MyVector(sympy_vector) or MyVector(sympy_matrix)
    
    ALGEBRAIC OPERATORS:
    --------------------
    v + w      # Vector addition
    v - w      # Vector subtraction  
    3 * v      # Scalar multiplication (returns vector)
    v * 3      # Scalar multiplication
    v.dot(w)   # Dot product (returns scalar)
    v @ w      # Dot product using @ operator
    v ** w     # Cross product using ** operator
    v | w      # Angle between vectors using | operator
    
    PROPERTIES:
    -----------
    • v.components    → (x, y, z) tuple
    • v.magnitude     → vector length
    • v.unit_vector   → normalized version
    • v.origin        → application point
    • v.endpoint      → terminal point
    • v.direction_cosines → (cosα, cosβ, cosγ)
    
    PHYSICS METHODS:
    ----------------
    • v.torque(reference_point) → torque vector τ = r × F
    • v.moment_about_axis(point, axis) → scalar moment about axis
    • v.projection(other) → vector projection
    
    EXAMPLES:
    ---------
    >>> v = MyVector(1, 2, 3)
    >>> w = MyVector(4, 5, 6)
    >>> 3 * v                     # MyVector(3, 6, 9)
    >>> v.dot(w)                  # 32
    >>> v ** w                    # MyVector(-3, 6, -3)
    >>> fuerza = MyVector([1,1,0], F=10, origin=(1,2,3))
    >>> torque = fuerza.torque()  # τ = r × F
    
    NOTES:
    ------
    • All operations support symbolic expressions via SymPy
    • Force vectors require F parameter and origin point
    • Coordinate system is shared globally across instances
    """
    
    def __init__(self, *args, origin=(0, 0, 0), F=None, dire=None):
        """
        Constructs a 3D vector with flexible input options.
        
        PARAMETERS:
        -----------
        *args : various
            Flexible input (see class docstring for options)
        origin : tuple/Point3D, default=(0,0,0)
            Vector application point for physics calculations
        F : number/symbol, optional
            Force magnitude (creates force vector when specified)
        dire : tuple/list, optional  
            Direction vector (normalized when F is provided)
        
        RAISES:
        -------
        ValueError : invalid argument combination
        TypeError : unsupported argument type
        
        EXAMPLES:
        ---------
        >>> MyVector(1, 2, 3)                     # From components
        >>> MyVector([4, 5, 6])                   # From list
        >>> MyVector(Point3D(0,0,0), Point3D(1,1,1)) # From points
        >>> MyVector([1,0,0], F=10, origin=(0,0,1))  # Force vector
        >>> MyVector(dire=(1,1,0), F=5)           # Force with direction
        """
         
        self.F = F
        self.dire = dire
        
        # Create coordinate system if it doesn't exist globally
        if not hasattr(MyVector, '_coord_sys'):
            MyVector._coord_sys = CoordSys3D('C')
        
        self.coord_sys = MyVector._coord_sys
        self.i = self.coord_sys.i
        self.j = self.coord_sys.j
        self.k = self.coord_sys.k
        
        # Set origin - CONVERTIR A SYMPY DIRECTAMENTE
        if isinstance(origin, Point3D):
            self._origin = Matrix([origin.x, origin.y, origin.z])
        else:
            origin = tuple(sympify(o) for o in origin)
            self._origin = Matrix(origin)
        
        # CASO ESPECIAL: Dos Point3D como argumentos
        if len(args) == 2 and all(isinstance(arg, Point3D) for arg in args):
            P1, P2 = args
            # El primer punto es el origen
            self._origin = Matrix([P1.x, P1.y, P1.z])
            
            # Calcular componentes como P2 - P1 - SIN FLOAT!
            dx = P2.x - P1.x
            dy = P2.y - P1.y
            dz = P2.z - P1.z
            
            self._components = Matrix([dx, dy, dz])
            self._parametric = False
            
            # Crear sympy vector
            self._sympy_vector = dx * self.i + dy * self.j + dz * self.k
            
            # Calcular endpoint
            self._endpoint = self._origin + self._components
            
            # Si hay F, escalar el vector para ajustar magnitud
            if F is not None:
                current_magnitude = self.magnitude
                if current_magnitude != 0:
                    scale_factor = F / current_magnitude
                    self._components = self._components * scale_factor
                    # Recrear sympy vector con componentes escaladas
                    dx, dy, dz = self._components
                    self._sympy_vector = dx * self.i + dy * self.j + dz * self.k
                    # Recalcular endpoint
                    self._endpoint = self._origin + self._components
            
            # Force properties
            self.is_force = F is not None
            self.force_magnitude = F
            self.force_components = self._components if self.is_force else None
            
            return
        
        # Si se proporciona dirección
        if dire is not None:
            if len(dire) != 3:
                raise ValueError("Direction vector must have exactly three components.")

            dire_vector = Matrix([sympify(d) for d in dire])
            dire_magnitude = sqrt(dire_vector.dot(dire_vector))

            if dire_magnitude == 0:
                raise ValueError("Direction vector cannot be zero")

            if F is not None:
                scale_factor = F / dire_magnitude
                self._components = dire_vector * scale_factor
                self._parametric = False
            else:
                t = symbols('t')
                self._components = dire_vector * t
                self._parametric = True
                 
                
        # Si no hay dirección, procesar args normalmente
        else:
            self._parametric = False
            
            # Process arguments
            if len(args) == 1:
                arg = args[0]
                if isinstance(arg, (list, tuple)):
                    # From list: [1, 2, 3]
                    self._components = Matrix([sympify(c) for c in arg])
                elif isinstance(arg, SympyVector):
                    # From sympy vector
                    self._components = Matrix([arg.dot(self.i), arg.dot(self.j), arg.dot(self.k)])
                elif isinstance(arg, Matrix):
                    # From sympy matrix
                    self._components = arg
                elif isinstance(arg, Point3D):
                    # From a single Point3D (tratar como vector posición desde origen)
                    self._origin = Matrix([0, 0, 0])
                    self._components = Matrix([arg.x, arg.y, arg.z])
                elif hasattr(arg, '__iter__'):
                    # From any iterable
                    self._components = Matrix([sympify(c) for c in arg])
                else:
                    raise TypeError("Invalid argument")
                    
            elif len(args) == 3:
                # Three separate numbers: MyVector(1, 2, 3)
                self._components = Matrix([sympify(c) for c in args])
                
            elif len(args) == 0:
                if dire is None and F is None:
                    self._components = Matrix([0, 0, 0])
                elif dire is not None:
                    dire_vector = Matrix(dire)
                    dire_magnitude = sqrt(dire_vector.dot(dire_vector))
                    if dire_magnitude == 0:
                        raise ValueError("Direction vector cannot be zero.")
                    if F is None:
                        self._components = dire_vector
                    else:
                        scale_factor = F / dire_magnitude
                        self._components = dire_vector * scale_factor
                else:
                    raise ValueError("Invalid initialization.")
        # Create sympy vector
        if not hasattr(self, '_sympy_vector'):
            dx, dy, dz = self._components
            self._sympy_vector = dx * self.i + dy * self.j + dz * self.k
        
        # Calculate endpoint
        self._endpoint = self._origin + self._components
        
        # Force properties
        self.is_force = F is not None
        self.force_magnitude = F
        self.force_components = self._components if self.is_force else None
        
        # Si F se proporcionó pero NO se usó dire, escalar para ajustar magnitud
        if self.is_force and not self._parametric and dire is None:
            current_magnitude = self.magnitude
            if current_magnitude != 0:
                scale_factor = F / current_magnitude
                self._components = self._components * scale_factor
                # Recreate sympy vector
                dx, dy, dz = self._components
                self._sympy_vector = dx * self.i + dy * self.j + dz * self.k
                # Recalculate endpoint
                self._endpoint = self._origin + self._components

    # Agregar propiedad para el punto final como Point3D
    def set(self,**kwargs):
        if len(kwargs)>0:
            self.i=real_subs(self.i,**kwargs)
            self.j=real_subs(self.j,**kwargs)
            self.k=real_subs(self.k,**kwargs)
        return self     
    @property
    def endpoint(self):
        """Returns the endpoint as a Point3D"""
        from sympy.geometry import Point3D
        x, y, z = self._endpoint
        return Point3D(x, y, z)
    
    @property
    def magnitude(self):
        """Returns the magnitude of the vector"""
        return sqrt(self._components.dot(self._components))
    
    @property
    def components(self):
        """Returns the vector components as a tuple"""
        return tuple(self._components)
    
    @property
    def direction_cosines(self):
        """Returns the direction cosines (α, β, γ)"""
        mag = self.magnitude
        if mag == 0:
            return (0, 0, 0)
        x, y, z = self._components
        return (x/mag, y/mag, z/mag)
    
    @property
    def unit_vector(self):
        """Returns the unit vector in the same direction"""
        mag = self.magnitude
        if mag == 0:
            return MyVector(0, 0, 0)
        x, y, z = self._components
        return MyVector(x/mag, y/mag, z/mag)

    def _clean_repr(self):
        """Clean representation without asterisks - SIMPLIFICADA"""
        x, y, z = self._components
        
        # Versión simple que evita comparaciones simbólicas
        def format_term(coef, unit):
            from sympy import simplify, sympify
            coef = simplify(sympify(coef))
            
            if coef == 0:
                return None
            elif coef == 1:
                return f"+{unit}"
            elif coef == -1:
                return f"-{unit}"
            else:
                return f"+({coef}){unit}"
        
        terms = []
        
        # i term
        term_i = format_term(x, 'i')
        if term_i:
            terms.append(term_i)
        
        # j term
        term_j = format_term(y, 'j')
        if term_j:
            terms.append(term_j)
        
        # k term
        term_k = format_term(z, 'k')
        if term_k:
            terms.append(term_k)
        
        if not terms:
            return "0"
        
        # Unir y limpiar
        result = ''.join(terms)
        
        # Quitar el + inicial si existe
        if result.startswith('+'):
            result = result[1:]
        
        # Reemplazar "+-" por "-"
        result = result.replace('+-', '-')
        
        return result
 
    def __repr__(self):
        """Official string representation"""
        return self._clean_repr()

    def __str__(self):
        """String representation"""
        return self.__repr__()
    
    @property
    def origin(self):
        """Returns origin as tuple"""
        return tuple(self._origin)
    
    @origin.setter
    def origin(self, new_origin):
        """Sets a new origin"""
        if isinstance(new_origin, (list, tuple)) and len(new_origin) == 3:
            self._origin = Matrix(new_origin)
            # Recalculate endpoint
            self._endpoint = self._origin + self._components
        elif isinstance(new_origin, Point3D):
            self._origin = Matrix([new_origin.x, new_origin.y, new_origin.z])
            # Recalculate endpoint
            self._endpoint = self._origin + self._components
        else:
            raise ValueError("Origin must be a tuple, list of 3 elements, or Point3D")
    
    @property
    def x(self):
        """X component"""
        return self._components[0]
    
    @property
    def y(self):
        """Y component"""
        return self._components[1]
    
    @property
    def z(self):
        """Z component"""
        return self._components[2]
    

    

    
    @property
    def matrix(self):
        """Components as column matrix"""
        return self._components
    
    @property
    def sympy_vector(self):
        """Equivalent sympy vector"""
        return self._sympy_vector
    

    @property 
    def module(self):
        """Alias for magnitude"""
        return self.magnitude
    
    @property
    def norm(self):
        """Alias for magnitude"""
        return self.magnitude
    
    @property
    def vunit(self):
        """Unit vector in same direction"""
        mag = self.magnitude
        if mag == 0:
            return MyVector(0, 0, 0, origin=self.origin)
        return MyVector([c/mag for c in self._components], origin=self.origin)
    
    def __add__(self, other):
        """Vector addition"""
        if isinstance(other, MyVector):
            new_components = self._components + other._components
            return MyVector(new_components, origin=self.origin)
        else:
            try:
                other_vec = MyVector(other)
                return self + other_vec
            except:
                raise TypeError(f"Cannot add MyVector with {type(other)}")
    
    def __radd__(self, other):
        """Right addition"""
        return self.__add__(other)
    
    def __sub__(self, other):
        """Vector subtraction"""
        if isinstance(other, MyVector):
            new_components = self._components - other._components
            return MyVector(new_components, origin=self.origin)
        else:
            try:
                other_vec = MyVector(other)
                return self - other_vec
            except:
                raise TypeError(f"Cannot subtract MyVector with {type(other)}")
    
    def __rsub__(self, other):
        """Right subtraction"""
        if isinstance(other, MyVector):
            return other - self
        else:
            try:
                other_vec = MyVector(other)
                return other_vec - self
            except:
                raise TypeError(f"Cannot subtract {type(other)} from MyVector")
    
    def __mul__(self, other):
        """Multiplication: vector * scalar → scaled vector"""
        from sympy import sympify
        
        # CASO 1: Si es otro MyVector -> ERROR (usa .dot() o @)
        if isinstance(other, MyVector):
            raise TypeError(
                "Use v1.dot(v2) for dot product or v1 @ v2 for matrix multiplication. "
                "For cross product use v1.cross(v2) or v1 ** v2."
            )
        
        # CASO 2: Si es un escalar (número, variable, expresión SymPy)
        try:
            # Convertir a expresión SymPy
            scalar = sympify(other)
            
            # Multiplicar componentes
            new_components = self._components * scalar
            
            # Crear nuevo vector
            result = MyVector(new_components, origin=self.origin)
            
            # Preservar propiedades de fuerza
            if hasattr(self, 'is_force') and self.is_force:
                result.is_force = True
                if hasattr(self, 'F') and self.F is not None:
                    result.F = self.F * scalar
                    result.force_magnitude = result.F
            
            return result
            
        except Exception as e:
            # Si no se pudo convertir a sympify, dar error claro
            raise TypeError(f"Cannot multiply MyVector by {type(other)}: {e}")
    
    def __rmul__(self, other):
        """Right multiplication - VERSIÓN CORREGIDA"""
        # La multiplicación escalar es conmutativa
        return self.__mul__(other)
    
    def __getitem__(self, index):
        """Component access by index"""
        if 0 <= index < 3:
            return self._components[index]
        elif index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError(f"Index {index} out of range. Vectors have 3 components")
    
    def dot(self, other):
        """Dot product"""
        if isinstance(other, MyVector):
            return self._components.dot(other._components)
        else:
            try:
                other_vec = MyVector(other)
                return self._components.dot(other_vec._components)
            except:
                raise TypeError(f"Cannot calculate dot product with {type(other)}")
    
    def cross(self, other):
        """Cross product"""
        if isinstance(other, MyVector):
            x1, y1, z1 = self._components
            x2, y2, z2 = other._components
            
            # Cross product: (y1*z2 - z1*y2, z1*x2 - x1*z2, x1*y2 - y1*x2)
            result = Matrix([
                y1*z2 - z1*y2,
                z1*x2 - x1*z2,
                x1*y2 - y1*x2
            ])
            
            return MyVector(result, origin=self.origin)
        else:
            try:
                other_vec = MyVector(other)
                return self.cross(other_vec)
            except:
                raise TypeError(f"Cannot calculate cross product with {type(other)}")

    def __matmul__(self, other):
        """Matrix multiplication operator @ for dot product"""
        return self.dot(other)

    def __pow__(self, other):
        """Cross product using ** operator"""
        return self.cross(other)
    
    def angle(self, other):
        """Angle between two vectors"""
        if isinstance(other, MyVector):
            dot_product = self.dot(other)
            mag_product = self.magnitude * other.magnitude
            if mag_product == 0:
                raise ValueError("Cannot calculate angle with null vector")
            return acos(dot_product / mag_product)
        else:
            try:
                other_vec = MyVector(other)
                return self.angle(other_vec)
            except:
                raise TypeError(f"Cannot calculate angle with {type(other)}")
    
    def __or__(self, other):
        """Angle using | operator"""
        return self.angle(other)
    
    def projection(self, other):
        """Projection of this vector onto another"""
        if isinstance(other, MyVector):
            dot_product = self.dot(other)
            other_mag_sq = other.magnitude**2
            if other_mag_sq == 0:
                return MyVector(0, 0, 0, origin=self.origin)
            factor = dot_product / other_mag_sq
            return MyVector([c * factor for c in other._components], origin=self.origin)
        else:
            try:
                other_vec = MyVector(other)
                return self.projection(other_vec)
            except:
                raise TypeError(f"Cannot calculate projection onto {type(other)}")
    
    def is_parallel(self, other, tolerance=1e-10):
        """Checks if two vectors are parallel"""
        if isinstance(other, MyVector):
            cross_product = self.cross(other)
            return abs(cross_product.magnitude) < tolerance
        else:
            try:
                other_vec = MyVector(other)
                return self.is_parallel(other_vec, tolerance)
            except:
                raise TypeError(f"Cannot check parallelism with {type(other)}")
    
    def is_orthogonal(self, other, tolerance=1e-10):
        """Checks if two vectors are orthogonal"""
        if isinstance(other, MyVector):
            dot_product = self.dot(other)
            return abs(dot_product) < tolerance
        else:
            try:
                other_vec = MyVector(other)
                return self.is_orthogonal(other_vec, tolerance)
            except:
                raise TypeError(f"Cannot check orthogonality with {type(other)}")
    
    def rotate(self, angle, axis='z'):
        """Rotates vector around an axis"""
        from sympy import cos, sin
        
        c, s = cos(angle), sin(angle)
        x, y, z = self._components
        
        if axis.lower() == 'x':
            # Rotation around x-axis
            new_components = Matrix([x, y*c - z*s, y*s + z*c])
        elif axis.lower() == 'y':
            # Rotation around y-axis
            new_components = Matrix([x*c + z*s, y, -x*s + z*c])
        elif axis.lower() == 'z':
            # Rotation around z-axis
            new_components = Matrix([x*c - y*s, x*s + y*c, z])
        else:
            raise ValueError("Axis must be 'x', 'y' or 'z'")
        
        return MyVector(new_components, origin=self.origin)
    
    def to_list(self):
        """Converts to list"""
        return self.components
    
    def to_tuple(self):
        """Converts to tuple"""
        return tuple(self._components)
    
    def to_matrix(self):
        """Converts to matrix"""
        return self._components.copy()
    
    def to_sympy_vector(self):
        """Converts to sympy vector"""
        return self._sympy_vector
        
    def s(self):
        display(Math(latex(self.components)))
        
    def display(self):
        """Displays complete vector information"""
        print(f"Vector: {self}")
        print(f"Origin: {self.origin}")
        print(f"Endpoint: {self.endpoint}")
        print(f"Components: {self.components}")
        print(f"Magnitude: {self.magnitude}")
        print(f"Unit vector: {self.vunit}")
        if self.is_force:
            print(f"Force magnitude: {self.force_magnitude}")
        return self
 
    
    def torque(self, reference_point=(0, 0, 0)):
        """
        Calculates torque VECTOR: τ = r × F
        
        Returns torque as a VECTOR
        """
        if not self.is_force:
            raise ValueError("Torque can only be calculated for force vectors")
        
        # Vector posición r
        rx = self.origin[0] - reference_point[0]
        ry = self.origin[1] - reference_point[1]
        rz = self.origin[2] - reference_point[2]
        
        r = MyVector([rx, ry, rz])
        
        # Torque vector τ = r × F
        τ_vector = r.cross(self)
        
        # Add properties
        τ_vector.is_torque = True
        τ_vector.torque_vector = τ_vector
        τ_vector.units = "N·m" if hasattr(self, 'units') else "force_units·m"
        τ_vector.reference_point = reference_point
        
        return τ_vector
    
    def moment(self, reference_point=(0, 0, 0), axis=None):
        """
        Calculates MOMENT MAGNITUDE (scalar) for equilibrium equations
        
        If axis is None: returns magnitude of torque vector
        If axis is given: returns moment about that specific axis
        
        Parameters:
        -----------
        reference_point: point about which to calculate moment
        axis: direction vector (for moment about an axis)
        
        Returns:
        --------
        float: moment magnitude
        """
        # First get torque vector
        τ = self.torque(reference_point)
        
        if axis is None:
            # Return magnitude of torque vector
            return τ.magnitude
        else:
            # Project torque onto axis to get moment about that axis
            if isinstance(axis, MyVector):
                axis_vector = axis.unit_vector
            else:
                axis_vector = MyVector(axis).unit_vector
            
            # Moment about axis = (τ · axis_unit) = projection
            moment_about_axis = τ.dot(axis_vector)
            return moment_about_axis
    
    def moment_about_axis(self, point_on_axis=(0, 0, 0), axis_direction=(0, 0, 1)):
        """
        Specific method for moment about an axis
        M = u · (r × F) where u is unit vector along axis
        """
        # Unit vector along axis
        u = MyVector(axis_direction).unit_vector
        
        # Vector from point on axis to force application point
        r = MyVector([
            self.origin[0] - point_on_axis[0],
            self.origin[1] - point_on_axis[1],
            self.origin[2] - point_on_axis[2]
        ])
        
        # Triple scalar product: u · (r × F)
        moment = u.dot(r.cross(self))
        return moment

        
# Utility functions
def zero_vector(origin=(0, 0, 0)):
    """Returns a zero vector"""
    return MyVector(0, 0, 0, origin=origin)

def vunit(axis='x', origin=(0, 0, 0)):
    """Returns a unit vector along specified axis"""
    if axis.lower() == 'x':
        return MyVector(1, 0, 0, origin=origin)
    elif axis.lower() == 'y':
        return MyVector(0, 1, 0, origin=origin)
    elif axis.lower() == 'z':
        return MyVector(0, 0, 1, origin=origin)
    else:
        raise ValueError("Axis must be 'x', 'y' or 'z'")

def from_sympy(sympy_vector, origin=(0, 0, 0)):
    """Creates MyVector from sympy vector"""
    return MyVector(sympy_vector, origin=origin)

def from_points(p1, p2):
    """Creates vector from point p1 to point p2"""
    if len(p1) != 3 or len(p2) != 3:
        raise ValueError("Points must have 3 coordinates")
    
    components = [p2[i] - p1[i] for i in range(3)]
    return MyVector(components, origin=p1)

def parallelogram_area(v1, v2):
    """Area of parallelogram formed by two vectors"""
    return v1.cross(v2).magnitude

def parallelepiped_volume(v1, v2, v3):
    """Volume of parallelepiped formed by three vectors"""
    return abs(v1.dot(v2.cross(v3)))
    
    
def solvevector(vector,*args):
    vec1=[vector.components[i] for i in range(3)]
    for data in args:
        vec1.append(data)
    return simplesolve(*vec1)    
        