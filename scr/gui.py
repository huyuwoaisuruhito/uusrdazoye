import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

from tkinter.scrolledtext import ScrolledText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import scr.file_io as fio
import scr.ddd_plot as dp
import scr.molecule as mol
import scr.toplevels as top

import scr.mm.computing as comp

from scr.mm.computing import DataTypeError
from scr.mm.atomtype import AtomTypeError

class Main_windows(tk.Tk):

    '''主窗口类'''

    def __init__(self):
        super().__init__()
        self.wm_title("测试")
        #self.geometry("1000x600")

        self.Computer = None
        self.Molecule = mol.Molecule()
        self.fio = fio.File_IO(self, self.Molecule)
        self.fio.input_gaussian_file('Gaussian_inp\\flat.gjf')#testing
        self.Molecule.auto_set_O()
        
        self.__menu = _Main_menu(self)
        self.__dd_frame = _DD_frame(self)
        self.ddd = _DDD_windows(self)
        self._log = Log(self)
        self.columnconfigure(0, weight=3)
        self.columnconfigure(1, weight=1)
        self.ddd.grid(row=0, column=0, sticky=tk.N+tk.S+tk.W+tk.E)
        self._log.grid(row=0, column=1, sticky=tk.N+tk.S+tk.W+tk.E)

        #self.Molecule.test()

        tk.mainloop()
    
    def open_bond_length_windows(self):
        if '_blw' in self.__dict__:
            self._blw.destroy()
        self._blw = top._Bond_windows(self)
    
    def open_bond_angle_windows(self):
        if '_baw' in self.__dict__:
            self._baw.destroy()
        self._baw = top._Bond_angle_windows(self)
    
    def open_dihedral_angle_windows(self):
        if '_daw' in self.__dict__:
            self._daw.destroy()
        self._daw = top._Dihedral_angle_windows(self)
    
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
    
    def warning(self, t):
        messagebox.showinfo('警告', t)


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
        editmenu.add_command(label="打开3D界面", command=self.__call_3d)
        editmenu.add_command(label="自动调整原点", command=self.__auto_set_O)
        editmenu.add_separator()
        editmenu.add_command(label="键长", command=self.parent.open_bond_length_windows)
        editmenu.add_command(label="键角", command=self.parent.open_bond_angle_windows)
        editmenu.add_command(label="二面角", command=self.parent.open_dihedral_angle_windows)

        #计算
        compmenu = tk.Menu(self.menubar, tearoff=0)
        compmenu.add_command(label="用MM优化分子构象", command=self.__com_MM)

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

    def __com_MM(self):
        self.parent.Computer = comp.Computing(self.parent, self.parent.Molecule)
        try:
            self.parent.Computer.run()
        except (AtomTypeError, DataTypeError) as error:
            self.parent.warning(error)

    def __auto_set_O(self):
        self.parent.Molecule.auto_set_O()
        self.parent.ddd.re_plot(self.parent.Molecule)

    def __call_3d(self):
        self.parent.open_3d_windows()

    def help_about(self):
        messagebox.showinfo('关于', '作者：文亦质 \n 1800011702 \n verion 0.1')

    def help_operate(self):
        messagebox.showinfo()



class _DDD_windows(tk.Frame):

    '''matplotlib_3d图像容器类'''

    def __init__(self, parent):
        super().__init__(parent)
        #self.wm_title("3D视图")
        #self.geometry("800x600")
        self._parent = parent

        self.plot = dp.DDD_plot()

        canvas = FigureCanvasTkAgg(self.plot.fig, self)
        canvas.draw()
        canvas_tkw = canvas.get_tk_widget()
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()

        toolbar.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=False)
        canvas_tkw.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.plot.init(self._parent.Molecule)
        # https://blog.csdn.net/qq_28485501/article/details/85329343 <- refer

        self._parent.bind('<Key-Up>', self.plot.change_view_pos)
        self._parent.bind('<Key-Down>', self.plot.change_view_pos)
        self._parent.bind('<Key-Left>', self.plot.change_view_pos)
        self._parent.bind('<Key-Right>', self.plot.change_view_pos)
        self._parent.bind('<MouseWheel>', self.plot.change_view_dist)
        self._parent.bind('<<Upadate3DView>>', self.view_call_back)

        self.plot.fig.canvas.mpl_connect('pick_event', self.__pick)
        self.plot.fig.canvas.mpl_connect('button_release_event', self.__button_release)

    def __pick(self, e):
        if e.mouseevent.button == 1:
            a = self.plot.heigh_light_atom(e.artist)
            if a>=0:
                self._parent.select(a)
            else:
                self._parent.clear_selection()
    
    def __button_release(self, e):
        if e.button == 3:
            self.plot.clear_high_light()
            self._parent.clear_selection()

    def re_plot(self, molecule):
        self.plot.re_plot(molecule)
    
    def view_call_back(self, event):
        self.re_plot(self._parent.Molecule)


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
        self._parent = parent
        
        self._parent.bind('<<UpadateLog>>', self.log_call_back)

    def log_call_back(self, event):
        self.__text['state'] = tk.NORMAL
        self.__text.insert(tk.END, '\n' + self._parent.Computer.message)
        self.__text.see(tk.END)
        self.__text['state'] = tk.DISABLED
        self.update()

