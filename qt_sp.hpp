
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

  /// We need to temporarly store this, since the computation takes a
  /// couple of seconds and would freeze the ui other wise. So we need
  /// to compute it using a signal (running it in the background
  /// thread), then we can directly access this value updating the UI
  double perc_thres;

  
signals:
  
  void parameter_changed();

  void geometry_changed();
  
  void projection_changed();
  
  void status_updated( QString st );

  void measurements_updated();
  
  void send_message( QString msg, int type = 1 ); 

  void set_status( int what, int state );					       


					
public slots:  

  double get_percolation_threshold();
  
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
