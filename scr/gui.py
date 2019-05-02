import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

from tkinter.scrolledtext import ScrolledText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import scr.file_io as fio
import scr.ddd_plot as dp
import scr.molecule as mol

Flag = True

class Main_windows(tk.Tk):

    '''主窗口类'''

    def __init__(self):
        super().__init__()
        self.wm_title("测试")
        self.geometry("800x600")

        self.__menu = _Main_menu(self)
        self.__dd_frame = _DD_frame(self)
        self.__log = Log(self)
        self.columnconfigure(0, weight=3)
        self.columnconfigure(1, weight=1)
        #self.__dd_frame.grid(row=0, column=1, sticky=tk.N+tk.S+tk.W)
        self.__log.grid(row=0, column=1, sticky=tk.N+tk.S+tk.E)

        self.Molecule = mol.Molecule()
        self.fio = fio.File_IO(self, self.Molecule)

        self.open_3d_windows()
        #self.Molecule.test()

        tk.mainloop()

    def open_3d_windows(self):
        self.ddd = _DDD_windows(self)
    
    def open_bond_length_windows(self):
        if '_blw' in self.__dict__:
            self._blw.destroy()
        self._blw = _Bond_length_windows(self)
    
    def open_bond_angle_windows(self):
        if '_baw' in self.__dict__:
            self._baw.destroy()
        self._baw = _Bond_angle_windows(self)
    
    def open_dihedral_angle_windows(self):
        if '_daw' in self.__dict__:
            self._daw.destroy()
        self._daw = _Dihedral_angle_windows(self)
    
    def select(self, a):
        if '_blw' in self.__dict__:
            self._blw.select(a)
        if '_baw' in self.__dict__:
            self._baw.select(a)
        if '_daw' in self.__dict__:
            self._daw.select(a)

    def clear_selection(self):
        if '_blw' in self.__dict__:
            self._blw.clear_selection()
        if '_baw' in self.__dict__:
            self._baw.clear_selection()
        if '_daw' in self.__dict__:
            self._daw.clear_selection()


class _Main_menu:

    '''菜单类'''
    
    def __init__(self, parent):

        '''初始化菜单'''

        self.menubar = tk.Menu(parent)
        self.parent = parent
        self.parent.config(menu=self.menubar)
        
        #文件
        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="打开", command=self.file_open)
        filemenu.add_command(label="新建", command=self.file_new)
        filemenu.add_command(label="保存", command=self.file_save)
        filemenu.add_separator()
        filemenu.add_command(label="退出", command=parent.quit)
        
        #编辑
        editmenu = tk.Menu(self.menubar, tearoff=0)
        editmenu.add_command(label="打开3D展示界面", command=self.call_3d)
        editmenu.add_separator()
        editmenu.add_command(label="键长", command=self.parent.open_bond_length_windows)
        editmenu.add_command(label="键角", command=self.parent.open_bond_angle_windows)
        editmenu.add_command(label="二面角", command=self.parent.open_dihedral_angle_windows)

        #计算
        compmenu = tk.Menu(self.menubar, tearoff=0)
        compmenu.add_command(label="用MM优化分子构象", command=self.call_3d)

        #帮助
        helpmenu = tk.Menu(self.menubar, tearoff=0)
        helpmenu.add_command(label="关于", command=self.help_about)
        helpmenu.add_command(label="使用方法", command=self.help_operate)
        
        self.menubar.add_cascade(label="文件", menu=filemenu)
        self.menubar.add_cascade(label="页面", menu=editmenu)
        self.menubar.add_cascade(label="计算", menu=compmenu)
        self.menubar.add_cascade(label="帮助", menu=helpmenu)

    def file_open(self):
        pass
    
    def file_new(self):
        pass
    
    def file_save(self):
        pass

    def _com_MM(self):
        pass

    def call_3d(self):
        self.parent.open_3d_windows()

    def help_about(self):
        messagebox.showinfo('关于', '作者：文亦质 \n 1800011702 \n verion 0.1')

    def help_operate(self):
        messagebox.showinfo()



class _Bond_length_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.wm_title("键长")
        self.resizable(width=False, height=False)
        self.Molecule = self.parent.Molecule.copy()
        self.__select = []
        
        self.a = tk.IntVar()
        self.b = tk.IntVar()
        self.l = tk.DoubleVar()
        self.l.trace("w", lambda name, index, mode: self.__change_text())

        f1 = ttk.Frame(self, padding=5); f1.grid(row=0, column=0)
        l1 = ttk.Label(f1, text='（不动的）第一个原子的序号：'); l1.grid(row=0, column=0)
        l2 = ttk.Label(f1, text='（可调的）第二个原子的序号：'); l2.grid(row=1, column=0)
        l3 = ttk.Label(f1, text='键长：'); l3.grid(row=2, column=0)

        e1 = ttk.Entry(f1, textvariable=self.a); e1.grid(row=0, column=1)
        e2 = ttk.Entry(f1, textvariable=self.b); e2.grid(row=1, column=1)
        e3 = ttk.Entry(f1, textvariable=self.l); e3.grid(row=2, column=1)

        f2 = ttk.Frame(self, padding=5); f2.grid(row=1, column=0)
        self.s1 = tk.Scale(f2, label='键长（Å）:', resolution=0.001, from_=0.001, to=4, orient=tk.HORIZONTAL, length=250, command=self.__change_scale)
        self.s1.grid(row=0, column=0)

        f3 = ttk.Frame(self, padding=5); f3.grid(row=2, column=0)
        b1 = ttk.Button(f3, text='确定', command=self.__commit); b1.grid(row=3, column=0)
        b2 = ttk.Button(f3, text='取消', command=self.__quit); b2.grid(row=3, column=1)

    def __change_scale(self, event):
        if self.l.get() != event:
            self.l.set(event)

    def __change_text(self):
        try:
            if self.s1.get() != self.l.get():
                self.s1.set(self.l.get())
                self.Molecule.modify_bond_length(self.a.get(), self.b.get(), self.l.get())
                self.parent.ddd.re_plot(self.Molecule)
        except RuntimeWarning:
            global Flag
            if Flag:
                Flag = False
                messagebox.showinfo('警告', 'RuntimeWarning，请检查输入')   
                def ___a():
                    global Flag
                    Flag = True
                self.after(50, ___a())

    def __commit(self):
        try:
                self.parent.Molecule.modify_bond_length(self.a.get(), self.b.get(), float(self.l.get()))
                self.parent.ddd.re_plot(self.parent.Molecule)
        except RuntimeWarning:
            global Flag
            if Flag:
                Flag = False
                messagebox.showinfo('警告', 'RuntimeWarning，请检查输入')   
                def ___a():
                    global Flag
                    Flag = True
                self.after(50, ___a())
        self.destroy()
    
    def __quit(self):
        self.parent.ddd.plot.clear_high_light()
        self.parent.ddd.re_plot(self.parent.Molecule)
        self.destroy()
        del self.parent._blw
    
    def __set_init(self):
        self.l.set(self.parent.Molecule.get_bond_length(self.a.get(), self.b.get()))
        self.s1.set(self.l.get())
    
    def select(self, i):
        self.__select.append(i)
        if len(self.__select) == 1:
            self.a.set(i)
        elif len(self.__select) == 2:
            self.b.set(i)
            self.__set_init()
        elif len(self.__select) >= 2:
            self.a.set(self.__select[-2])
            self.b.set(self.__select[-1])
            self.__set_init()
    
    def clear_selection(self):
        self.__select = []
        self.a.set(0)
        self.b.set(0)
        self.l.set(0)

class _Bond_angle_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.wm_title("键角")
        self.resizable(width=False, height=False)
        self.Molecule = self.parent.Molecule.copy()
        self.__select = []

        self.o = tk.IntVar()
        self.a = tk.IntVar()
        self.b = tk.IntVar()
        self.angle = tk.DoubleVar()
        self.angle.trace("w", lambda name, index, mode: self.__change_text())

        f1 = ttk.Frame(self, padding=5); f1.grid(row=0, column=0)
        l1 = ttk.Label(f1, text='（不动的）第一个原子的序号：'); l1.grid(row=0, column=0)
        l2 = ttk.Label(f1, text='（不动的）顶点原子的序号：'); l2.grid(row=1, column=0)
        l3 = ttk.Label(f1, text='（可调的）第三个原子的序号：'); l3.grid(row=2, column=0)
        l4 = ttk.Label(f1, text='键角：'); l4.grid(row=3, column=0)

        e1 = ttk.Entry(f1, textvariable=self.a); e1.grid(row=0, column=1)
        e2 = ttk.Entry(f1, textvariable=self.o); e2.grid(row=1, column=1)
        e3 = ttk.Entry(f1, textvariable=self.b); e3.grid(row=2, column=1)
        e4 = ttk.Entry(f1, textvariable=self.angle); e4.grid(row=3, column=1)

        f2 = ttk.Frame(self, padding=5); f2.grid(row=1, column=0)
        self.s1 = tk.Scale(f2, label='键角:', resolution=0.1, from_=0, to=360, orient=tk.HORIZONTAL, length=250, command=self.__change_scale)
        self.s1.grid(row=0, column=0)

        f3 = ttk.Frame(self, padding=5); f3.grid(row=2, column=0)
        b1 = ttk.Button(f3, text='确定', command=self.__commit); b1.grid(row=4, column=0)
        b2 = ttk.Button(f3, text='取消', command=self.__quit); b2.grid(row=4, column=1)

    def __change_scale(self, event):
        if self.angle.get() != event:
            self.angle.set(event)

    def __change_text(self):
        try:
            if self.s1.get() != self.angle.get():
                self.s1.set(self.angle.get())
                self.Molecule.modify_bond_angle(self.a.get(), self.o.get(), self.b.get(), self.angle.get()/180*np.pi)
                self.parent.ddd.re_plot(self.Molecule)  
        except RuntimeWarning:
            global Flag
            if Flag:
                Flag = False
                messagebox.showinfo('警告', 'RuntimeWarning，请检查输入')   
                def ___a():
                    global Flag
                    Flag = True
                self.after(50, ___a())

    def __commit(self):
        try:
            self.parent.Molecule.modify_bond_angle(self.a.get(), self.o.get(), self.b.get(), self.angle.get()/180*np.pi)
            self.parent.ddd.re_plot(self.parent.Molecule)
            self.destroy()
        except RuntimeWarning:
            global Flag
            if Flag:
                Flag = False
                messagebox.showinfo('警告', 'RuntimeWarning，请检查输入')   
                def ___a():
                    global Flag
                    Flag = True
                self.after(50, ___a())
    
    def __quit(self):
        self.parent.ddd.plot.clear_high_light()
        self.parent.ddd.re_plot(self.parent.Molecule)
        self.destroy()
        del self.parent._baw

    def __set_init(self):
        angle = self.parent.Molecule.get_bond_angle(self.a.get(), self.o.get(), self.b.get())/np.pi*180
        self.angle.set(self.parent.Molecule.get_bond_angle(self.a.get(), self.o.get(), self.b.get())/np.pi*180)
        self.s1.set(self.angle.get())

    def select(self, i):
        self.__select.append(i)
        if len(self.__select) == 1:
            self.a.set(i)
        elif len(self.__select) == 2:
            self.o.set(i)
        elif len(self.__select) == 3:
            self.b.set(i)
            self.__set_init()
        elif len(self.__select) >= 3:
            self.a.set(self.__select[-3])
            self.o.set(self.__select[-2])
            self.b.set(self.__select[-1])
            self.__set_init()
    
    def clear_selection(self):
        self.__select = []
        self.a.set(0)
        self.o.set(0)
        self.b.set(0)
        self.angle.set(0)


class _Dihedral_angle_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.wm_title("二面角")
        self.resizable(width=False, height=False)
        self.Molecule = self.parent.Molecule.copy()
        self.__select = []

        self.a = tk.IntVar()
        self.b = tk.IntVar()
        self.c = tk.IntVar()
        self.d = tk.IntVar()
        self.angle = tk.DoubleVar()
        self.angle.trace("w", lambda name, index, mode: self.__change_text())

        f1 = ttk.Frame(self, padding=5); f1.grid(row=0, column=0)
        l1 = ttk.Label(f1, text='（不动的）第一个原子的序号：'); l1.grid(row=0, column=0)
        l2 = ttk.Label(f1, text='（不动的）第二个原子的序号：'); l2.grid(row=1, column=0)
        l3 = ttk.Label(f1, text='（不动的）第三个原子的序号：'); l3.grid(row=2, column=0)
        l4 = ttk.Label(f1, text='（可调的）第四个原子的序号：'); l4.grid(row=3, column=0)
        l5 = ttk.Label(f1, text='二面角：'); l5.grid(row=4, column=0)

        e1 = ttk.Entry(f1, textvariable=self.a); e1.grid(row=0, column=1)
        e2 = ttk.Entry(f1, textvariable=self.b); e2.grid(row=1, column=1)
        e3 = ttk.Entry(f1, textvariable=self.c); e3.grid(row=2, column=1)
        e4 = ttk.Entry(f1, textvariable=self.d); e4.grid(row=3, column=1)
        e5 = ttk.Entry(f1, textvariable=self.angle); e5.grid(row=4, column=1)

        f2 = ttk.Frame(self, padding=5); f2.grid(row=1, column=0)
        self.s1 = tk.Scale(f2, label='二面角:', resolution=0.1, from_=-180, to=180, orient=tk.HORIZONTAL, length=250, command=self.__change_scale)
        self.s1.grid(row=0, column=0)

        f3 = ttk.Frame(self, padding=5); f3.grid(row=2, column=0)
        b1 = ttk.Button(f3, text='确定', command=self.__commit); b1.grid(row=5, column=0)
        b2 = ttk.Button(f3, text='取消', command=self.__quit); b2.grid(row=5, column=1)

    def __change_scale(self, event):
        if self.angle.get() != event:
            self.angle.set(event)

    def __change_text(self):
        try:
            if self.s1.get() != self.angle.get():
                self.s1.set(self.angle.get())
                self.Molecule.modify_dihedral_angle(self.a.get(), self.b.get(), self.c.get(), self.d.get(), self.angle.get()/180*np.pi)
                self.parent.ddd.re_plot(self.Molecule)
        except RuntimeWarning:
            global Flag
            if Flag:
                Flag = False
                messagebox.showinfo('警告', 'RuntimeWarning，请检查输入')   
                def ___a():
                    global Flag
                    Flag = True
                self.after(50, ___a())

    def __commit(self):
        try:
            self.parent.Molecule.modify_dihedral_angle(self.a.get(), self.b.get(), self.c.get(), self.d.get(), self.angle.get()/180*np.pi)
            self.parent.ddd.re_plot(self.parent.Molecule)
            self.destroy()
        except RuntimeWarning:
            global Flag
            if Flag:
                Flag = False
                messagebox.showinfo('警告', 'RuntimeWarning，请检查输入')   
                def ___a():
                    global Flag
                    Flag = True
                self.after(50, ___a())
    
    def __quit(self):
        self.parent.ddd.plot.clear_high_light()
        self.parent.ddd.re_plot(self.parent.Molecule)
        self.destroy()
        del self.parent._daw

    def __set_init(self):
        self.angle.set(self.parent.Molecule.get_dihedral_angle(self.a.get(), self.b.get(), self.c.get(), self.d.get())/np.pi*180)
        self.s1.set(self.angle.get())

    def select(self, i):
        self.__select.append(i)
        if len(self.__select) == 1:
            self.a.set(i)
        elif len(self.__select) == 2:
            self.b.set(i)
        elif len(self.__select) == 3:
            self.c.set(i)
        elif len(self.__select) == 4:
            self.d.set(i)
            self.__set_init()
        elif len(self.__select) >= 4:
            self.a.set(self.__select[-4])
            self.b.set(self.__select[-3])
            self.c.set(self.__select[-2])
            self.d.set(self.__select[-1])
            self.__set_init()
    
    def clear_selection(self):
        self.__select = []
        self.a.set(0)
        self.b.set(0)
        self.c.set(0)
        self.d.set(0)
        self.angle.set(0)



class _DDD_windows(tk.Toplevel):

    '''matplotlib_3d图像容器类'''

    def __init__(self, parent):
        super().__init__(parent)
        self.wm_title("3D视图")
        self.geometry("800x600")
        self.parent = parent

        self.plot = dp.DDD_plot()
        self.parent.fio.input_gaussian_file('Gaussian_inp\\HCHO.gjf')#testing

        canvas = FigureCanvasTkAgg(self.plot.fig, self)
        canvas.draw()
        canvas_tkw = canvas.get_tk_widget()
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()

        toolbar.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=False)
        canvas_tkw.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.plot.init(self.parent.Molecule)
        # https://blog.csdn.net/qq_28485501/article/details/85329343 <- refer

        self.bind('<Key-Up>', self.plot.change_view_pos)
        self.bind('<Key-Down>', self.plot.change_view_pos)
        self.bind('<Key-Left>', self.plot.change_view_pos)
        self.bind('<Key-Right>', self.plot.change_view_pos)
        self.bind('<MouseWheel>', self.plot.change_view_dist)

        self.plot.fig.canvas.mpl_connect('pick_event', self.__pick)
        self.plot.fig.canvas.mpl_connect('button_release_event', self.__button_release)

    def __pick(self, e):
        if e.mouseevent.button == 1:
            a = self.plot.heigh_light_atom(e.artist)
            if a>=0:
                self.parent.select(a)
            else:
                self.parent.clear_selection()
    
    def __button_release(self, e):
        if e.button == 3:
            self.plot.clear_high_light()
            self.parent.clear_selection()

    def re_plot(self, molecule):
        self.plot.re_plot(molecule)



class _DD_frame(tk.Frame):

    '''2D绘图类'''

    def __init__(self, parent):
        super().__init__(parent)

    def height_light_atom(self, a):
        pass

    def height_light_bond(self, a, b):
        pass



class Log(tk.Frame):

    '''计算过程的Log画面'''

    def __init__(self, parent):
        super().__init__(parent)

        self.__label = tk.Label(self, text='Log')
        self.__text = tk.scrolledtext.ScrolledText(self, width=40, height=40, state = tk.DISABLED)
        self.__label.grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        self.__text.grid(row=1, column=0, sticky=tk.N+tk.S+tk.E+tk.W)

    def log_call_back(self, event):
        self.__text['state'] = tk.NORMAL
        self.__text.insert(tk.END, '\n' + event)#event不能提供信息，交互有待实现
        self.__text['state'] = tk.DISABLED