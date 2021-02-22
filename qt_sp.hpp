/* Projection tool - compute planar projection of triply periodic
 * minimal surfaces 
 * Copyright (C) 2020 Tobias Hain
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see https://www.gnu.org/licenses.
 */



#ifndef QT_SP_H
#define QT_SP_H

#include <QTimer>
#include <QObject>
#include <QApplication>
#include <QThread>
#include <QFont>
#include <QFontMetrics>
#include <QImage>
#include <QPen>
#include <QPainter>

#include <iostream>
#include <iomanip>
#include <sstream>

#include "surface_projection.hpp"


/** \brief A wrapper class around surface_projection to use with a QT interface
 *
 */
class sp_qt : public QObject, public surface_projection {

  Q_OBJECT;
  
public:
  sp_qt ();
  ~sp_qt();


  sp_qt( const sp_qt& ) = default;
  
private:


  
signals:
  
  void parameter_changed();

  void geometry_changed();
  
  void projection_changed();
  
  void status_updated( QString st );

  void measurements_updated( QString what );
  
  void send_message( QString msg, int type = 1 ); 

  void set_status( int what, int state );					       


					
public slots:  

  uchar* save_png_legend( unsigned char* img,
			  unsigned int width,
			  unsigned int height,
			  std::string textcolor,
			  std::string backcolor,
			  std::string fn,
			  std::vector<std::string> parameters = std::vector<std::string> (0, "") );
  
  void update_geometry_();  
  void compute_projection();

  void update_measurements( QString what );
  
  void change_surface_type( int ind );

  void change_uc_scale_ab( double ab );
  void change_uc_scale_c( double c );

  void set_slice_dim_to_uc();
  
  void change_vol_prop( double val );
  void change_lvl_set( double val );
  void change_xy_points( int val );
  void change_z_points( int val );
  void change_hkl( int h, int k, int l );
  void change_slice_thickness( double val );
  void change_slice_height( double val );
  void change_slice_width( double val );  
  void change_slice_position( double val );  
  void change_membranes( std::vector<double> val );

  void do_something();

  void copy_parameters( sp_qt *source );

  void save_grid( QString fn );
  void save_topological_network( int id, QString fn_prefix );
  void save_surface_points( int id, QString fn_prefix );

  void change_channel_color( int id, int val );
  
};


#endif
