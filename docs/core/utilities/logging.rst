Logging
=======

The Acts logging facility supports several severity levels which allow you
to control the amount of information displayed at run-time. Logger objects
can easily be created using the :func:`Acts::getDefaultLogger` function which
should be sufficient to get you started. In case you need more customized
debug output, you can make use of the output decorators defined in
:any:`Acts::Logging` or even write your own implementation of
:class:`Acts::Logging::OutputDecorator`. In order to add debug messages to your program, you should use the provided macros for the different severity levels:

.. code-block::

   ACTS_VERBOSE(...);
   ACTS_DEBUG(...);
   ACTS_INFO(...);
   ACTS_WARNING(...);
   ACTS_ERROR(...);
   ACTS_FATAL(...);

All of these macros require that a function ``logger()`` returning a
:class:`Acts::Logger` object is available in the scope in which the macros are
used. Inside classes containing an :class:`Acts::Logger` object as member
variable, this could be achieved by providing a private class method called
``logger()`` (for an example see e.g.
:func:`Acts::CylinderVolumeBuilder::logger`). Inside free functions or member
methods with local logger objects, the same effect can be achieved by using the
macro ``ACTS_LOCAL_LOGGER(...)`` which is provided for your convenience.

Code example illustrating the usage:

.. code-block::

   #include <fstream>
   #include <memory>
   
   #include "Acts/Utilities/Logger.hpp"
   
   void myFunction() {
     // open the logfile
     std::ofstream logfile("log.txt");
     // setup a logger instance for >= INFO messages, streaming into the log file
     // make sure you do NOT call the variable 'logger'
     std::unique_ptr<const Acts::Logger> myLogger
         = Acts::getDefaultLogger("MyLogger", Acts::Logging::INFO, &logfile);
     // make sure the Acts debug macros can work with your logger
     ACTS_LOCAL_LOGGER(myLogger);
     ACTS_VERBOSE("This message will not appear in the logfile.");
     ACTS_INFO("But this one will: Hello World!");
     // do not forget to close the logfile
     logfile.close();
   }

In case you are using Acts in another framework which comes with its own
logging facility (e.g. Gaudi) you can pipe the logging output from Acts
tools and algorithms to your framework's logging system by supplying different
implementations of:

- :class:`Acts::Logging::OutputFilterPolicy` (for mapping logging levels)
- :class:`Acts::Logging::OutputPrintPolicy` (for passing the Acts output
    to your internal logging system)

Since Acts makes extensive use of :func:`Acts::getDefaultLogger()` to provide
sufficient information for debugging, you would need to provide a modified
implementation of this function (using your output filter and printing policies)
to also pipe this output to your framework. Changing the implementation of an
already defined function is a non-trivial task in C++. We recommend the
following approach using the possibility to inject custom code by pre-loading
shared libraries with ``LD_PRELOAD``. You need to provide an appropriate
implementation for a function of the following signature into a separate source
file and compile it in a shared library

.. code-block::

   namespace Acts {
   std::unique_ptr<const Logger> getDefaultLogger(const std::string&,
                                                  const Logging::Level&,
                                                  std::ostream*);
   }

Then you can run your executable, which uses Acts tools and algorithms, in
the following way (tested under Unix)

.. code-block:: console

   LD_PRELOAD=<YOUR_SHARED_LIBRARY> path/to/your/exectuable

For an example have a look at CustomDefaultLogger.cpp which you can use as
follows:

.. code-block:: console

   cd <ACTS/INSTALL/DIRECTORY>
   source bin/setup.sh
   LD_PRELOAD=lib/libActsCustomLogger.so bin/Examples/ActsGenericDetector
