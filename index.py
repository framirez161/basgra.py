from flask import Flask, render_template, render_template_string, send_file, request, redirect, Response, url_for, flash

from BASGRA import mainBasgra
from xlrd import open_workbook  
from matplotlib import pyplot 

import os

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, static_url_path = "/tpm", static_folder = "tpm")

@app.route('/')
def home():
    return render_template('map.html')
   
@app.route('/result', methods=['GET', 'POST'])
def result():    
    path_1 = ''
    path_1 = APP_ROOT + "\\tpm\\out.png"
    path_1 = path_1.replace('\\','/')
   
    try:
        os.remove(path_1)
    except:
        pass    
    
    
    projectpath = request.form['projectFilepath']
    split_ =  projectpath.split(' ,')[0]
    split__ = split_.split(',')[0]
    split___ = split__.replace('.', '')
    area = float(split___)

    """
    Open and read a parameters set in Excel file
    """
    
    path_1 = ''
    path_1 = APP_ROOT + "\\static\\parameters_test.xlsx"
    book = open_workbook(path_1.replace('\\','/'))
    sheet = book.sheet_by_index(0)   
    param = []
    for i in range (sheet.nrows):
        param.append([])
        for j in range (1,5):
            cells = sheet.row_slice(rowx=i,start_colx=j,end_colx=j+1)
            for cell in cells:
                if i > 0:
                    param[i].append(float(cell.value)) 
                else:
                    param[i].append(cell.value)
   
    """
    Open and read a weather set in Excel file
    """
    path_1 = ''
    path_1 = APP_ROOT + "\\static\\weather_test.xlsx"
    book = open_workbook(path_1.replace('\\','/'))
    sheet = book.sheet_by_index(0)   
    weather = []
    for i in range (sheet.nrows):
        weather.append([])
        for j in range (10):
            cells = sheet.row_slice(rowx=i,start_colx=j,end_colx=j+1)
            for cell in cells:
                weather[i].append(cell.value)

    #print(weather[1][1])
        
    #mainBasgra(param,weather,DAYS_HARVEST,NDAYS,NOUT,y)
    
      #EVERY 35 DAYS
    DAYS_HARVEST = [[2017,144],[2017,179],[2017,214]]
       
    asign = [0,1,2,3] # asing parameters
    x = []
    starday = weather[1][2]
      
    for j in range(sheet.nrows-1):
        x.append(starday)
        starday = starday + 1
            
    y = []
    for i in range(sheet.nrows-1):
        y.append([])

    dm_ = []        
    for i in range(len(asign)):
        out = mainBasgra(param,weather,asign[i],DAYS_HARVEST)    
        dm = out[28] # grass dry matter
        dm_.append([])
        dm_[i] = dm
        
        for j in range(len(dm)):
            y[j].append((dm[j]/1000) * area) 

    pyplot.xlabel("Dia del año (2017)")
    pyplot.ylabel("Producción de MS de pasto (Kg/m2)")
    #pyplot.title("GRASS REGROWTH IN THREE LEVELS  (Area: " + str(projectpath)+')')
    labels = ['Priori','Alto','Medio','Bajo']
    for i in range(len(labels)):
        pyplot.plot(x,[pt[i] for pt in y],label = labels[i])
    pyplot.legend()
    path_1 = ''
    path_1 = APP_ROOT + "\\tpm\\out.png"
    #dd= APP_ROOT + "\\tpm\\out.png"
    path_1 = path_1.replace('\\','/')
    pyplot.savefig(path_1, dpi= 200)
    pyplot.show(block=False)
    pyplot.close()

    return render_template('result.html',dm=dm_[0])

# 
#----------------------------------------------------------------------

if __name__ == '__main__':
    app.run(debug=True)