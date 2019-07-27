'''
Created on 16/03/2011

@author: ispmarin
'''
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Polygon, Ellipse, Circle
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
import numpy as np
import math
import sys

import aux_functions
import Crack
import LineDouble
import CircleInhom
import ErrorAnalysis


class plotter(object):
    '''
    classdocs
    '''
    def __init__(self, sx1, sy1, sx2, sy2, sstep):
        '''
        Constructor
        '''
        self.x1 = sx1
        self.x2 = sx2
        self.y1 = sy1
        self.y2 = sy2
        self.step = sstep
        self.high_head_lines = 100
        self.low_head_lines = 40
        self.streamlines = 50
        #print( 'Plotter: x_min= %s,  x_max= %s,  y_min= %s, y_max= %s, num_steps= %s' % (self.x1, self.x2, self.y1, self.y2, self.step))
    
    def plot_potential(self, crack_list, u_flow, element_name, save_to_file):
        '''
        Plot the potential
        '''
        print('plotting...')
        dim_x = np.linspace(self.x1, self.x2, self.step)
        dim_y = np.linspace(self.y1, self.y2, self.step)
        X, Y = np.meshgrid(dim_x, dim_y)
        
        len_x = len(dim_x)
        len_y = len(dim_y)
        
        
        a = np.zeros([len_x, len_y], dtype=complex)
              
        for i, y in enumerate(dim_y):
            for j, x in enumerate(dim_x):
                a[i][j] = aux_functions.final_potential(complex(x, y), crack_list, u_flow)
                
        
        plt.rcParams['contour.negative_linestyle'] = 'solid'
        
        formatter = FuncFormatter(lambda x,pos: ("%.2e"%x).replace(".",","))
        fig0 = plt.figure(0)

        fig0.suptitle("Potencial")
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        cs = plt.contourf(X, Y, (a.real), self.low_head_lines)
        cs2 = plt.contour(X, Y, (a.real), self.low_head_lines, linewidths=0.6, colors='black')
        plt.clabel(cs2, inline=1,fontsize=8)
        plt.colorbar(cs,format=formatter)
        ax = plt.gca()
        self.add_polygon(crack_list,ax)
        if save_to_file == 'yes':
            save_name1 = str.lower(element_name) + '_pot_'+str(self.low_head_lines) + '.pdf'
            plt.savefig(save_name1,  facecolor='w', edgecolor='w',
                        format='pdf',bbox_inches='tight')

        fig1 = plt.figure(1)
        fig1.suptitle("Linhas de Corrente")
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.contour(X, Y, a.imag, self.streamlines, linewidths=0.6)
        ax = plt.gca()
        self.add_polygon(crack_list,ax)
        if save_to_file == 'yes':
            save_name1 = str.lower(element_name) + '_stream_' +str(self.streamlines)+ '.pdf'
            plt.savefig(save_name1,  facecolor='w', edgecolor='w',
                        format='pdf',bbox_inches='tight')
        
        fig2 = plt.figure(2)
        fig2.suptitle("Potencial")
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        cs = plt.contourf(X, Y, (a.real), self.high_head_lines)
        cs2 = plt.contour(X, Y, (a.real), self.high_head_lines, linewidths=0.6, colors='black')
        ax = plt.gca()
        self.add_polygon(crack_list,ax)
        if save_to_file == 'yes':
            save_name1 = str.lower(element_name) + '_pot_' +str(self.high_head_lines)+ '.pdf'
            plt.savefig(save_name1,  facecolor='w', edgecolor='w',
                        format='pdf',bbox_inches='tight')
        
        #plt.show()
        
        print("done.")
        
    def plot_head(self, crack_list, u_flow, aq, element_name, kind_solver, save_to_file):
        '''
        Plot the potential
        '''
        
        
        
        print('plotting...')
        dim_x = np.linspace(self.x1, self.x2, self.step)
        dim_y = np.linspace(self.y1, self.y2, self.step)
        X, Y = np.meshgrid(dim_x, dim_y)

        
        len_x = len(dim_x)
        len_y = len(dim_y)
      
        a = np.zeros([len_x, len_y], dtype=complex)
        b = np.zeros([len_x, len_y], dtype=complex)
              
        for i, y in enumerate(dim_y):
            for j, x in enumerate(dim_x):
                a[i][j] = aux_functions.head_value(complex(x, y), crack_list, u_flow,aq)
                b[i][j] = aux_functions.final_potential(complex(x, y), crack_list, u_flow)
                
        plt.rcParams['contour.negative_linestyle'] = 'solid'
        
        formatter = FuncFormatter(lambda x,pos: ("%.2e"%x).replace(".",","))

        fig0 = plt.figure(0)
        #Head, values with lines
        fig0.suptitle("Carga (m)")
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        cs = plt.contourf(X, Y, (a.real), self.low_head_lines)
        cs2 = plt.contour(X, Y, (a.real), self.low_head_lines, linewidths=0.6, colors='black')
        plt.clabel(cs2, inline=1,fontsize=8)
        plt.colorbar(cs,format=formatter)
        ax = plt.gca()
        self.add_polygon(crack_list,ax)
        if save_to_file == 'yes':
            save_name1 = str.lower(element_name) + '_head_' + str(self.low_head_lines) + '_'+kind_solver + '.pdf'
            plt.savefig(save_name1,  facecolor='w', edgecolor='w',
                        format='pdf',bbox_inches='tight')
        
        plt.clf()

        fig1 = plt.figure(1)
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        #Streamlines, solid black
        fig1.suptitle("Linhas de Corrente")
        #cs3 = plt.contourf(X, Y, a.imag, 800)
        plt.contour(X, Y, b.imag, self.streamlines, linewidths=0.6)
        ax = plt.gca()
        self.add_polygon(crack_list,ax)
        if save_to_file == 'yes':
            save_name1 = str.lower(element_name) + '_stream_'+str(self.streamlines) +'_'+kind_solver+ '.pdf'
            plt.savefig(save_name1,  facecolor='w', edgecolor='w',
                        format='pdf',bbox_inches='tight')
        plt.clf()
        
        
        fig2 = plt.figure(2)
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        #Head, only lines
        fig2.suptitle("Carga (m)")
        cs4 = plt.contourf(X, Y, (a.real), 150)
        plt.contour(X, Y, (a.real), self.high_head_lines, linewidths=0.6, colors='black')
        plt.colorbar(cs4,format=formatter)
        ax = plt.gca()
        self.add_polygon(crack_list,ax)
        if save_to_file == 'yes':
            save_name1 = str.lower(element_name) + '_head_' + str(self.high_head_lines) +'_'+kind_solver + '.pdf'
            plt.savefig(save_name1,  facecolor='w', edgecolor='w',
                        format='pdf',bbox_inches='tight')
        plt.clf()
        
        fig3 = plt.figure(3)
        fig3.set_facecolor('#0099cc')
        fig3.set_edgecolor('#0099cc')
        plt.xlabel('X')
        plt.ylabel('Y')
        #Head, only lines
        #fig3.suptitle("Carga (m)")
        #cs4 = plt.contourf(X, Y, (a.real), 50)
        plt.contour(X, Y, (a.real), self.low_head_lines, linewidths=0.6, colors='white')
        plt.contour(X, Y, b.imag, self.streamlines, linewidths=0.6,linestyles='--', colors='white',alpha=0.6)
        #plt.colorbar(cs4,format=formatter)
        ax = plt.gca()
        ax.patch.set_facecolor('#0099cc')
        self.add_polygon(crack_list,ax)
        if save_to_file == 'yes':
            save_name1 = str.lower(element_name) + '_head_plain_' + str(self.high_head_lines) +'_'+kind_solver + '.pdf'
            plt.savefig(save_name1,  facecolor='#0099cc', edgecolor='#0099cc',
                        format='pdf',bbox_inches='tight')   
  
        #plt.show()
        plt.clf()
        
        print("done.")
          
        
    def plot_potential_matplotlib(self, crack_list, line_double_list, u_flow, aq,
                                  absolute_tolerance, relative_tolerance, kind_compare, 
                                  name_sol_1,name_sol_2, kind_solver, save_to_file):
        '''
        Plot the potential
        '''
        formatter = FuncFormatter(lambda x,pos: ("%.2e"%x).replace(".",","))
        
        fontP = FontProperties()
        fontP.set_size('small')

        print( 'plotting...')
        dim_x = np.linspace(self.x1, self.x2, self.step)
        dim_y = np.linspace(self.y1, self.y2, self.step)
        X, Y = np.meshgrid(dim_x, dim_y)
        #Z = X +1j*Y
        
        len_x = len(dim_x)
        len_y = len(dim_y)
        
        #print len_x, len_y
        
        a = np.zeros([len_x, len_y], dtype=complex)
        b = np.zeros([len_x, len_y], dtype=complex)
        
              
        for i, y in enumerate(dim_y):
            for j, x in enumerate(dim_x):
                if kind_compare == 'head':
                    a[i][j] = aux_functions.head_value(complex(x, y), crack_list, u_flow, aq)
                    b[i][j] = aux_functions.head_value(complex(x, y), line_double_list, u_flow,aq)
                elif kind_compare == 'potential':
                    a[i][j] = aux_functions.final_potential(complex(x, y), crack_list, u_flow)
                    b[i][j] = aux_functions.final_potential(complex(x, y), line_double_list, u_flow)
                else:
                    sys.exit('ERROR Please choose between head or potential in comparison')
                    
        
        absolute_error = ErrorAnalysis.compare_results_absolute(a, b)
        relative_error = ErrorAnalysis.compare_results_relative(a, b)
        
        file_err_name = 'err_' + kind_solver + '.dat' 
        f = open(str(file_err_name),'w')
        
        max_err = 'rel ' + str(np.max(relative_error)) + '\n'
        f.write(str(max_err))
        max_mean_err = 'mean relative error' + str(np.sum(relative_error)/(relative_error.shape[0] * relative_error.shape[1])) + '\n', 
        f.write(str(max_mean_err))
        
        max_err = 'abs ' + str(np.max(absolute_error)) + '\n'
        f.write(str(max_err))
        max_mean_err = 'mean absolute error' + str(np.sum(absolute_error)/(absolute_error.shape[0] * absolute_error.shape[1])) + '\n', 
        f.write(str(max_mean_err))
#        fig0 = plt.figure(0)
        #Head, first element list
        #plt.subplot(211)
#        if kind_compare == 'head':
#            fig0.suptitle('Carga (m)')
#        else:
#            fig0.suptitle('Potencial')
#            
#        cs = plt.contourf(X, Y, a.real, 800)
#        plt.contour(X, Y, a.real, self.high_head_lines, linewidths=0.6, colors='black')
#        plt.colorbar(cs,format="%.2f")
#        plt.xlabel('X (m)')
#        plt.ylabel('Y (m)')
#        ax = plt.gca()
#        self.add_polygon(crack_list,ax)
#        if save_to_file == 'yes':
#            save_name1 = str.lower(name_sol_1) +''+str(self.high_head_lines)+'_'+ kind_solver + '.pdf'
#            plt.savefig(save_name1,  facecolor='w', edgecolor='w',
#                        format='pdf')
#        
#        fig1 = plt.figure(1)
#        #plt.subplot(212)
#        if kind_compare == 'head':
#            fig1.suptitle('Carga (m)')
#        else:
#            fig1.suptitle('Potencial')
#        cs3 = plt.contourf(X, Y, b.real, 800)
#        plt.contour(X, Y, b.real, self.high_head_lines, linewidths=0.6, colors='black')
#        plt.colorbar(cs3,format="%.2f")
#        plt.xlabel('X (m)')
#        plt.ylabel('Y (m)')
#        ax = plt.gca()
#        self.add_polygon(line_double_list,ax)
#        if save_to_file == 'yes':
#            save_name2 = str.lower(name_sol_2) +'_'+ str(self.high_head_lines) +'_'+ kind_solver + '.pdf'
#            plt.savefig(save_name2,  facecolor='w', edgecolor='w',
#                        format='pdf')
#       

        if np.max(absolute_error) > absolute_tolerance:
            fig2 =plt.figure(2)
            fig2.suptitle("Erro Absoluto")
            cs3 = plt.contourf(X, Y, absolute_error.real, 300)
            plt.contour(X, Y, absolute_error.real, self.high_head_lines, linewidths=0.6, colors='black')
            plt.colorbar(cs3,format=formatter)
            if save_to_file == 'yes':
                save_name_abs = str.lower(name_sol_1) + '_' + str.lower(name_sol_2) + '_'+ 'abs_'+ kind_solver+ '.pdf'
                plt.savefig(save_name_abs,  facecolor='w', edgecolor='w',
                        format='pdf', bbox_inches='tight')
            plt.clf()
            
        if np.max(relative_error) > relative_tolerance:
            fig3 = plt.figure(3)
            fig3.suptitle("Erro Relativo")
            cs4 = plt.contourf(X, Y, relative_error.real, 300)
            plt.contour(X, Y, relative_error.real, self.high_head_lines, linewidths=0.6, colors='black')
            plt.colorbar(cs4,format=formatter)
            if save_to_file == 'yes':
                save_name_rel = str.lower(name_sol_1) + '_' + str.lower(name_sol_2) +  '_rel_'+ kind_solver+ '.pdf'
                plt.savefig(save_name_rel,  facecolor='w', edgecolor='w',
                        format='pdf', bbox_inches='tight')

            plt.clf()
        
        fig4 = plt.figure(4)
        fig4.suptitle("Comparativo")
        ax2 = fig4.add_subplot(111)#plt.gca()
        plot1 = ax2.contour(X, Y, a.real, self.low_head_lines, linewidths=0.8, colors = 'g', linestyles='dotted')
        plot2 = ax2.contour(X, Y, b.real, self.low_head_lines, linewidths=0.6, colors = 'b', linestyles='dashed')
        lines1 = plot1.collections[0]
        lines2 = plot2.collections[0]
        ax2.legend([lines1,lines2],(name_sol_1, name_sol_2),'upper right', shadow=True,prop=fontP)
        

        self.add_polygon(crack_list,ax2)
        if save_to_file == 'yes':
            save_name = str.lower(name_sol_1) + '_' + str.lower(name_sol_2) +  '_' + kind_compare + '_' + kind_solver + '.pdf'
            plt.savefig(save_name,  facecolor='w', edgecolor='w',
                        format='pdf', bbox_inches='tight')
        
        #plt.show()
        
        print('done.')
        
        

    
    def plot_to_file(self, crack_list, u_flow):
        '''
        Plot the potential
        '''
        print( 'plotting...')
        x = np.linspace(self.x1, self.x2, self.step)
        y = np.linspace(self.y1, self.y2, self.step)

        
        f = open('potential.out', 'w')
        fi = open('streamline.out', 'w')
       
        for i in range(len(x)):
            for j in range(len(y)):
                stt = '  %8f      ' % (aux_functions.final_potential(complex(x[j], y[i]), crack_list, u_flow).real)
                stt2 = '  %8f      ' % (aux_functions.final_potential(complex(x[j], y[i]), crack_list, u_flow).imag)
                f.write(stt)
                fi.write(stt2)
            stt = '\n'
            f.write(stt)
            fi.write(stt)
        
        f.close()
        fi.close()
        print( 'done.')



    def plot_coeffs(self, filename,n):
        data = np.loadtxt(filename)
        x = np.zeros(1)
        y1= np.zeros(1)
        y2= np.zeros(1)

        for i in xrange(n-1):
            print 'i', i, 'i+1+n', i+1+n, 'n',n
            
            for line in data:
                x = np.append(x,line[0])
                y1 = np.append(y1, line[i + 1])
                y2 = np.append(y2, line[i + n  ])
            plt.figure(i)
            plt.plot(x,y1, x,y2)
            stro = 'coefficient %d' %i
            plt.title(stro)
            x = np.zeros(1)
            y1 = np.zeros(1)
            y2 = np.zeros(1)
        
        plt.show()
        
        
    def add_polygon(self,crack_list,ax):

        elem_id = 0
        same_elem_id = []
        local_vertices = []
        patches = []
        
        if type(crack_list[0]) == LineDouble.LineDouble:
        
            for element in crack_list:
                if element.id == elem_id:
                    local_vertices.append([element.z1.real, element.z1.imag])
                else:
                    same_elem_id.append(local_vertices)
                    local_vertices = []
                    local_vertices.append([element.z1.real, element.z1.imag])
                    elem_id += 1
            same_elem_id.append(local_vertices)
            
            for vertices_list in same_elem_id:
                polygon = Polygon(vertices_list, facecolor='g',
                edgecolor='white', linewidth=3, alpha=0.1)
                patches.append(polygon)
            
        elif type(crack_list[0]) == Crack.Crack:
            for element in crack_list:
                center = (element.z1 + element.z2)/2
                
                patches.append(Ellipse((center.real,center.imag), element.L, element.aperture, angle=180*element.angle/math.pi, facecolor='black',
                edgecolor='w', linewidth=3, alpha=0.1 ))
        elif type(crack_list[0]) == CircleInhom.CircleInhom or type(crack_list[0]) == CircleInhom.ExactSolution:
            for element in crack_list:
                center = element.center
                patches.append(Circle((center.real, center.imag), element.R, fc='g',
                edgecolor='w', linewidth=3, alpha=0.1))
        #p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.5)
        else:
            print("Drawing polygon not implemented!")
            
        p = PatchCollection(patches)
        p.set_alpha(0.3)
        p.set_facecolor('black')
        p.set_edgecolor('white')
        p.set_linewidth(3)
        
        ax.set_aspect('equal')
        ax.add_collection(p)

    def check_head_line(self,general_element_list, u_flow, pos_ini, pos_end):
        
        n_pts = 200
        dim_x = np.linspace(0, 1, n_pts)
        
        x = np.zeros(n_pts)
        y = [np.zeros(n_pts)  for _ in xrange(len(general_element_list))]
        
        x2x1 = pos_end.real - pos_ini.real
        y2y1 = pos_end.imag - pos_ini.imag
    
        for j, element_list in enumerate(general_element_list):
            for i,t in enumerate(dim_x):
                z = complex(pos_ini.real + t*x2x1, pos_ini.imag + t*y2y1) 
                x[i] = t
                y[j][i] = aux_functions.head_value(z, element_list, u_flow)
        
        fig = plt.figure(1)
        ax2 = fig.add_subplot(111)
        #for j in enumerate(general_element_list):
        plt.plot(x,y[0],'k--',x,y[1],'r:',x,y[2],'g',x,y[3],'b-.')
        plt.xlabel('Linha, normalizada')
        locs,labels = plt.yticks()
        plt.yticks(locs, map(lambda x: "%.2f" % x, locs))
        plt.ylabel('Carga')
        ax2.legend(('n=1','n=2', 'n=4','n=8'),'upper right', shadow=True)
        
        save_name = 'line_double_head_check' + '.pdf'
        plt.savefig(save_name,  facecolor='w', edgecolor='w',
                        format='pdf', bbox_inches='tight')
        
        #plt.show()
        