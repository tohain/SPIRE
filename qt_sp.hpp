
#ifndef QT_SP_H
#define QT_SP_H

#include <QTimer>
#include <QObject>
#include <QApplication>
#include <QThread>

#include "surface_projection.hpp"

class sp_qt : public QObject, public surface_projection {

  Q_OBJECT;
  
public:
  sp_qt (double &p, std::string &s);
  ~sp_qt();


  sp_qt( const sp_qt& ) = default;
  
private:
  
signals:


  void parameter_changed();

  void geometry_changed();
  
  void projection_changed();
  
  void status_updated( QString st );

  void measurements_updated();
  
  void send_message( QString msg, int type = 1 ); 

  void set_status( int what, int state );					       
						
public slots:

  void update_geometry_();  
  void compute_projection();

  void update_measurements();
  
  void change_surface_type( int ind );
  void change_ntucs( int ind );
  void change_uc_size_a( double a );
  void change_uc_size_c( double a );  
  void change_vol_prop( double val );
  void change_xy_points( int val );
  void change_z_points( int val );
  void change_hkl( int h, int k, int l );
  void change_slice_width( double val );
  void change_slice_position( double val );  
  void change_membranes( std::vector<double> val );

  void do_something();

  void copy_parameters( sp_qt *source );

  void save_grid( QString fn );
  void save_topological_network( int id, QString fn_prefix );
  void save_surface_points( int id, QString fn_prefix );
  
  
};


#endif
