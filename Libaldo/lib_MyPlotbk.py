import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from lib_MyEq import *
from lib_MyEqEq import *
from matplotlib.axis import Axis
from lib_tools import *
import mplcursors
from sympy import *
import warnings

def vecspace(x1, x2, x3):
    tt = np.linspace(float(x1), float(x2), x3)
    return tt

def dataplotXY(obj, x1, x2, x3=200):
    obj = obj2str(obj)
    vecx = vecspace(x1, x2, x3)
    xx, yy = dataplot(vecx, obj)
    xx1, xx2 = min(xx), max(xx)
    vecx = getinterx(obj, var=x)
    
    if len(vecx) > 0:
        ee = obj2MyEq(obj)
        xmin, xmax = min(vecx), max(vecx)
        if xmin < xx1:
            y1 = float(ee(x=xmin))
            xx = [xmin] + xx
            yy = [y1] + yy
        if xmax > xx2:
            y2 = float(ee(x=xmax))
            xx = xx + [xmax]
            yy = yy + [y2]
        
    if len(xx) < 100:
        return dataplotXY(obj, xx1, xx2, x3=200)
    else:    
        return xx, yy, xx1, xx2

def dataplot(vecx, sfunc):
    X = []
    Y = []
    ee = MyEq(parse_expr(sfunc), 'ee', var=x, show=False)
    for i in vecx:
        kres = ee(x=i)
        if not 'I' in str(kres):
            Y.append(ee(x=i))
            X.append(i)
    return X, Y

def obj2str(obj):
    if type(obj) == str:
        return obj
    elif type(obj) == MyEq:
        return str(obj.ksym)
    elif type(obj) == MyEqEq:
        return str(simplify(expand(obj.L - obj.R)))
    else:
        return str(obj)

def obj2func(obj):
    if type(obj) == str:
        return parse_expr(obj)
    elif type(obj) == MyEq:
        return obj.ksym 
    elif type(obj) == MyEqEq:
        return simplify(expand(obj.L - obj.R)) 
    else:
        return obj

def obj2MyEq(obj, var=x):
    obj = obj2func(obj)
    ee = MyEq(obj, 'ee', var=var, show=False)
    return ee

def getinterx(obj, var=x):
    ee = obj2MyEq(obj, var=var)
    vx = ee.roots(show=False)
    if 'I' in str(vx) and ee.degree() == 3:
        vec3 = ee.coef_list()
        vx = list(solve3(*vec3))
    vx2 = []
    for i in vx:
        if not 'I' in str(i):
            vx2.append(i)
    return vx2
 
class MyPlotC:
    def __init__(self, *args, x_min=None, x_max=None, y_min=None, y_max=None, samples=1000, tam=1.0):
        functions = []
        vecf = []
        
        for data in args:
            if type(data) == MyEq:
                vecf.append(data)
                functions.append(str(data.ksym))
            elif type(data) == MyEqEq:
                vecf.append(MyEq(data.e1.ksym - data.e2.ksym, 'ee', show=False))
                functions.append(str(data.e1.ksym - data.e2.ksym))
            else:
                vecf.append(MyEq(data, 'ee', show=False))
                functions.append(str(data))
                
        self.functions_str = functions
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.vecx = []
        self.vecy = []
        
        # 🔥 APLICAR FACTOR DE TAMAÑO
        self.size_factor = tam
        self.samples = int(samples * tam)  # Muestras proporcionales
        
        self.x = symbols('x')
        self.functions = [sympify(f) for f in functions]
        self.lambdas = [self._safe_lambdify(f) for f in self.functions]
        
        # Calcular rangos automáticos si no se especifican
        if self.x_min is None or self.x_max is None:
            self.x_min, self.x_max = self._calculate_x_range()
        
        self.x_vals = np.linspace(self.x_min, self.x_max, self.samples)
        self.y_vals = [lam(self.x_vals) for lam in self.lambdas]
        
        # 🔥 TAMAÑO DE FIGURA PROPORCIONAL
        base_width, base_height = int(12*tam), int(6*tam)
        figsize = base_width, base_height
        
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self._plot()

    def _calculate_x_range(self):
        """Calcula rango X automáticamente basado en las funciones"""
        # Rango por defecto amplio
        x_min, x_max = -10, 10
        
        try:
            # Buscar puntos interesantes (raíces, puntos críticos)
            all_points = []
            for f in self.functions:
                # Raíces
                try:
                    roots = solve(f, self.x)
                    for r in roots:
                        if r.is_real:
                            all_points.append(float(r))
                except:
                    pass
                
                # Puntos críticos
                try:
                    df = diff(f, self.x)
                    criticals = solve(df, self.x)
                    for c in criticals:
                        if c.is_real:
                            all_points.append(float(c))
                except:
                    pass
            
            if all_points:
                x_min = min(all_points) - 2
                x_max = max(all_points) + 2
                # Asegurar que el rango no sea demasiado pequeño
                if x_max - x_min < 4:
                    x_min, x_max = -5, 5
        except:
            pass
        
        return x_min, x_max

    def _safe_lambdify(self, f):
        f_lamb = lambdify(self.x, f, modules=["numpy"])
        def wrapper(x_array):
            try:
                y = f_lamb(x_array)
                y = np.asarray(y, dtype=np.complex128)
                y.real[np.abs(y.imag) > 1e-8] = np.nan
                return y.real
            except:
                return np.full_like(x_array, np.nan, dtype=float)
        return wrapper

    def _calculate_y_range(self):
        """Calcula rango Y automáticamente"""
        if self.y_min is not None and self.y_max is not None:
            return self.y_min, self.y_max
            
        all_y = []
        for y_vals in self.y_vals:
            valid_y = y_vals[~np.isnan(y_vals)]
            if len(valid_y) > 0:
                all_y.extend(valid_y)
        
        if len(all_y) == 0:
            return -10, 10
            
        y_min, y_max = np.min(all_y), np.max(all_y)
        padding = (y_max - y_min) * 0.1  # 10% de padding
        
        return y_min - padding, y_max + padding

    def _find_intersections(self):
        intersections = []
        for i in range(len(self.functions)):
            for j in range(i + 1, len(self.functions)):
                eq = Eq(self.functions[i], self.functions[j])
                sol = solve(eq, self.x)
                for s in sol:
                    if s.is_real:
                        try:
                            x_val = float(s.evalf())
                            if self.x_min <= x_val <= self.x_max:
                                y1 = self.functions[i].subs(self.x, x_val).evalf()
                                y2 = self.functions[j].subs(self.x, x_val).evalf()
                                if not y1.is_real or not y2.is_real:
                                    continue
                                y_val = float(y1)
                                intersections.append((x_val, y_val))
                                self.vecx.append(x_val)
                                self.vecy.append(y_val) 
                        except:
                            continue
        return intersections
 
    def _find_critical_points(self):
        critical_points = []
        for f in self.functions:
            df = diff(f, self.x)
            sol = solve(df, self.x)
            for s in sol:
                if s.is_real:
                    try:
                        s_val = float(s.evalf())
                        if self.x_min <= s_val <= self.x_max:
                            y_val = float(f.subs(self.x, s).evalf())
                            critical_points.append((s_val, y_val))
                            self.vecx.append(s_val)
                            self.vecy.append(y_val)
                    except:
                        continue
        return critical_points

    def _find_inflexion_points(self):
        inflection_points = []
        delta = (self.x_max - self.x_min) / self.samples
        for f in self.functions:
            d2f = diff(f, self.x, 2)
            sol = solve(d2f, self.x)
            for s in sol:
                if s.is_real:
                    try:
                        x0 = float(s.evalf())
                        if self.x_min <= x0 <= self.x_max:
                            left = d2f.subs(self.x, x0 - delta)
                            right = d2f.subs(self.x, x0 + delta)
                            if left * right < 0:
                                y_val = float(f.subs(self.x, x0).evalf())
                                inflection_points.append((x0, y_val))
                                self.vecx.append(x0)
                                self.vecy.append(y_val)
                    except:
                        continue
        return inflection_points

    def _find_axis_intersections(self):
        points = []
        for f in self.functions:
            # Eje Y
            try:
                y0 = float(f.subs(self.x, 0).evalf())
                if not np.isnan(y0) and np.isreal(y0):
                    points.append((0, y0))
                    self.vecy.append(y0)
            except:
                pass

            # Eje X
            sol = solve(f, self.x)
            for s in sol:
                if s.is_real:
                    try:
                        x0 = float(s.evalf())
                        if self.x_min <= x0 <= self.x_max:
                            points.append((x0, 0))
                            self.vecx.append(x0)
                    except:
                        continue
        return points    

    def _plot(self):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.functions)))
        criticals = []
        inflections = []
        intersections = []
        axis_inters = []
        
        # 1. GRAFICAR FUNCIONES (esto sí mantuve)
        for i, y in enumerate(self.y_vals):
            label = fr"$f_{{{i+1}}}(x) = {latex(self.functions[i])}$"
            self.ax.plot(self.x_vals, y, label=label, color=colors[i])

        # 2. 🔥 ESTO ES LO QUE QUITÉ - CALCULAR PUNTOS ESPECIALES
        try:
            criticals = self._find_critical_points()
        except:
            pass    

        try:
            inflections = self._find_inflexion_points()
        except:
            pass

        try:
            intersections = self._find_intersections()
        except:
            pass

        try:
            axis_inters = self._find_axis_intersections()
        except:
            pass        

        # 3. 🔥 ESTO TAMBIÉN QUITÉ - MOSTRAR PUNTOS ESPECIALES
        all_special = (
            [(x, y, "Máximo/Mínimo") for x, y in criticals] +
            [(x, y, "Inflexión") for x, y in inflections] + 
            [(x, y, "Intersección de curvas") for x, y in intersections] +
            [(x, y, "Eje X/Y") for x, y in axis_inters]
        )

        if all_special:
            xs, ys, labels = zip(*all_special)
            sizes = [20 if label == "Eje X/Y" else 40 for label in labels]
            scatter = self.ax.scatter(xs, ys, color='red', s=sizes, zorder=5)
            
            # Tooltips al hover
            try:
                cursor = mplcursors.cursor(scatter, hover=True)
                @cursor.connect("add")
                def on_hover(sel):
                    x, y = sel.target
                    idx = sel.index
                    tag = labels[idx]
                    sel.annotation.set(text=f"{tag}\nx = {x:.3f}\ny = {y:.3f}", backgroundcolor='white')
            except:
                pass

        # 4. CONFIGURACIÓN FINAL (esto sí mantuve)
        self.ax.axhline(0, color='black', lw=1.5)
        self.ax.axvline(0, color='black', lw=1.5)
        self.ax.grid(True)
        self.ax.legend()
        plt.show()
        
def MyPlot(*args, xmin=None, xmax=None, ymin=None, ymax=None, samples=1000,tam=1):
    '''
    MyPlot(f1, f2, Q1, expr)
    f1, f2 = MyEq
    Q1 = MyEqEq  
    expr = x*4+5...
    
    mejoras en Jup Lab
    ####################################################
    from IPython.display import display, HTML
    display(HTML("""
    <style>
        .container { width: 98% !important; }
        .output_result { max-width: 100% !important; }
        .jp-Cell-outputArea { overflow: visible !important; }
        #notebook-container { width: 99% !important; margin: 0 auto !important; }
    </style>
    """))

    import matplotlib.pyplot as plt
    plt.rcParams['figure.figsize'] = (24, 6)
    plt.rcParams['figure.dpi'] = 60
    ####################################################
    '''
 
 
    x1=-1
    x2=1
    samples=samples
    cc=0
    vdata=[]
    for data in args:
        if not Is_Number(data):
            vdata.append(data)
        elif cc==0:
            x1=data
            cc+=1
        elif cc==1:
            x2=data
            cc+=1
        else:
            samples=data
    
    # 🔥 Aplicar factor de tamaño
    if hasattr(MyPlotC, '_original_init'):
        # Restaurar init original si ya fue parcheado
        MyPlotC.__init__ = MyPlotC._original_init
    
    # Parchear temporalmente la clase para pasar el tamaño
    original_init = MyPlotC.__init__
    
    def new_init(self, *args, **kwargs):
        # Ajustar tamaño antes de llamar al init original
        if 'size' in kwargs:
            size_factor = kwargs.pop('size')
            # Ajustar samples proporcionalmente
            if 'samples' in kwargs:
                kwargs['samples'] = int(kwargs['samples'] * size_factor)
        original_init(self, *args, **kwargs)
    
    MyPlotC.__init__ = new_init
    MyPlotC._original_init = original_init
       
    PP = MyPlotC(*vdata, 
                 x_min=xmin, 
                 x_max=xmax, 
                 y_min=ymin,
                 y_max=ymax,
                 samples=int(samples * tam),  # 🔥 Ajustar samples
                 tam=tam)  # 🔥 Pasar factor de tamaño