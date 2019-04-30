import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

from tkinter.scrolledtext import ScrolledText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import file_io as fio
import ddd_plot as dp
import molecule as mol

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

    def open_3d_windows(self):
        self.__ddd = _DDD_windows(self)
    
    def open_bond_length_windows(self):
        _Bond_length_windows(self)
    
    def open_bond_angle_windows(self):
        _Bond_angle_windows(self)
    
    def open_dihedral_angle_windows(self):
        _Dihedral_angle_windows(self)




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
        editmenu = tk.Menu(self.menubar, tearoff=0)
        editmenu.add_command(label="计算MM", command=self.call_3d)

        #帮助
        helpmenu = tk.Menu(self.menubar, tearoff=0)
        helpmenu.add_command(label="关于", command=self.help_about)
        helpmenu.add_command(label="使用方法", command=self.help_operate)
        
        self.menubar.add_cascade(label="文件", menu=filemenu)
        self.menubar.add_cascade(label="页面", menu=editmenu)
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
        super().__init__()
        self.wm_title("键长")
        self.geometry("300x200")
        self.resizable(width=False, height=False)
        self.initial_windows()

    def initial_windows(self):
        pass

class _Bond_angle_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__()
        self.wm_title("键角")
        self.geometry("300x200")
        self.resizable(width=False, height=False)
        self.initial_windows()

    def initial_windows(self):
        pass


class _Dihedral_angle_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__()
        self.wm_title("二面角")
        self.geometry("300x200")
        self.resizable(width=False, height=False)
        self.initial_windows()

    def initial_windows(self):
        pass
        


class _DDD_windows(tk.Toplevel):

    '''matplotlib_3d图像容器类'''

    def __init__(self, parent):
        super().__init__(parent)
        self.wm_title("3D视图")
        self.geometry("800x600")
        self.parent = parent

        plot=dp.DDD_plot()
        self.parent.fio.Gaussian_file('Gaussian_inp\\CHCl3.gjf')#testing

        canvas = FigureCanvasTkAgg(plot.fig, self)
        canvas.draw()
        canvas_tkw = canvas.get_tk_widget()
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()

        toolbar.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=False)
        canvas_tkw.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        plot.init(self.parent.fio.atoms, self.parent.fio.bonding)
        # https://blog.csdn.net/qq_28485501/article/details/85329343 <- refer



class _DD_frame(tk.Frame):

    '''2D绘图'''

    def __init__(self, parent):
        super().__init__(parent)
        



class Log(tk.Frame):

    '''计算过程的Log画面'''

    def __init__(self, parent):
        super().__init__(parent)

        self.__label = tk.Label(self, text='Log')
        self.__text = tk.scrolledtext.ScrolledText(self, width=40, height=40, state = tk.DISABLED)
        self.__label.grid(row=0, column=0, in_=self, sticky=tk.N+tk.S+tk.E+tk.W)
        self.__text.grid(row=1, column=0, in_=self, sticky=tk.N+tk.S+tk.E+tk.W)

    def log_call_back(self, event, **kw):
        self.__text['state'] = tk.NORMAL
        self.__text.insert(tk.END, '\n' + event)
        self.__text['state'] = tk.DISABLED


if __name__ == "__main__":
    m = Main_windows()
    tk.mainloop()