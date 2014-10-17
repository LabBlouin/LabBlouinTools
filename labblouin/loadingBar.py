import Tkinter

class loadingBar:
        
    # Create Progress Bar
    def __init__(self, title, width=300, height=25):
        self.obj = Tkinter.Tk()
        self.obj.resizable(False, False)
        self.obj.title(title)
        self.canvas = Tkinter.Canvas(self.obj, width=width, height=height)
        self.canvas.grid()
        self.width = width
        self.height = height

    # Open Progress Bar
    def open(self):
        self.obj.deiconify()
        self.obj.focus_set()
        #self.obj.update()

    # Close Progress Bar
    def close(self):
        self.obj.withdraw()

    # Update Progress Bar
    def update(self, ratio):
        self.canvas.delete(Tkinter.ALL)
        self.canvas.create_rectangle(0, 0, self.width * ratio, \
                                       self.height, fill='blue')
        self.obj.update()
        self.obj.focus_set()
        
    # Set Title
    def settitle(self, title):
        self.obj.title(title)    