# Function to run on different operational systems

# This script provides helper functions to help
# detect the OS to run the correct code or save the outputs:
f.get_os = function(){
  sysinf = Sys.info()
  if (!is.null(sysinf)){
    os = sysinf['sysname']
    if (os == 'Darwin')
      os = "osx"
  } else { ## mystery machine
    os = .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os = "osx"
    if (grepl("linux-gnu", R.version$os))
      os = "linux"
  }
  tolower(os)
}

# get the current time and date to save the outputs:
f.get_date = function(){
  format(Sys.time(), "%y%m%d_%H%M")}