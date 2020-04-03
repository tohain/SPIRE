
#ifndef QT_SP_H
#define QT_SP_H


#include <QObject>
#include <QApplication>

#include "surface_projection.hpp"

class sp_qt : public QObject, public surface_projection {

  Q_OBJECT;
  
public:
  sp_qt (double &p, std::string &s);
  ~sp_qt();


private:


signals:


  void parameter_changed();

  void geometry_changed();
  
  void projection_changed();
  
  void status_updated( QString st );
			  
public slots:

  void update_geometry_();  
  void compute_projection();

  
  void change_surface_type( int ind );
  void change_ntucs( int ind );
  void change_uc_size( double a );
  void change_vol_prop( double val );
  void change_xy_points( int val );
  void change_z_points( int val );
  void change_h( int val );
  void change_k( int val );
  void change_l( int val );
  void change_slice_width( double val );
  void change_slice_position( double val );  
  void change_membranes( std::vector<double> val );
  
};


#endif
