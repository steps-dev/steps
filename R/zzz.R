# Check for OSGeo4W.bat and gdal_polygonize.py
.onLoad <-
  function(libname, pkgname) {
    op <-
      base::options()
    op.dlmpr <-
      base::list(
        cmd.path = {
          if (.Platform$OS.type == "windows") {
            base::Sys.which("OSGeo4W.bat")
          } else {
            'python'
          }
        },
        py.path = {
          if (.Platform$OS.type == "windows") {
            'gdal_polygonize'
            } else {
            base::Sys.which('gdal_polygonize.py')
          }
        }
      )
    toset <-
      op.dlmpr   %>%
      base::names(.) %>%
      magrittr::is_in(
        x     = .,
        table = base::names(op)
      ) %>%
      magrittr::not(.);
    if (base::any(toset)) {
      base::options(op.dlmpr[toset]);
    };
    base::invisible();
  }
