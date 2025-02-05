"""
MODULE FOR SET-UP, LOGGING AND DATE-OVERRIDE
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,platform,datetime,logging,builtins,time,multiprocessing,inspect

from importlib.metadata import version
from importlib import import_module


def set_date(date_toggle, date_override = datetime.date(2024,7,31)):
    if date_toggle == 0:
        today = datetime.date.today()
    elif date_toggle == 1:
        today = date_override

    actual_today = datetime.date.today()
    date_str = today.strftime('%Y-%m-%d')
    curr_year = int(str(today)[:4])
    overall_start_time = time.time()
    formatted_start_time = datetime.datetime.now().strftime("%H:%M:%S")

    return today, actual_today, date_str, curr_year, overall_start_time, formatted_start_time


def current_function_name():
    return inspect.currentframe().f_back.f_code.co_name


def import_and_log(*args):
    """
    Import a package or a component from a package and log it along 
    with its version immediately.
    """
    try:
        
        # Importing based on the number of arguments:
        # if 1 argument, 'import <arg_1>' (whole package)
        if len(args) == 1:
            package_name = args[0]
            package = import_module(package_name)
            imported_packages[package_name] = package
            setattr(builtins, package_name, package)  # globally available
            msg = f"Package '{package_name}' imported."
            
        # if 2 arguments, 'from <arg_1> import <arg_2>' (component from a package)
        elif len(args) == 2:
            package_name, component_name = args
            package = import_module(package_name)
            component = getattr(package, component_name)
            imported_packages[f"{package_name}.{component_name}"] = component
            setattr(builtins, component_name, component)  # globally available
            msg = f"Component '{component_name}' from package '{package_name}' imported."
        else:
            raise ValueError("Invalid number of arguments passed to import_and_log().")
        
        log_imports_and_versions(msg)  # Log the import and its version
        
    except ImportError as e:
        logging.error(f"Failed to import {args}. Error: {e}")
    except AttributeError as e:
        logging.error(f"Component '{args[1]}' not found in package '{args[0]}'. Error: {e}")


def log_imports_and_versions(message):
    """
    Append import messages and version information of the most
    recently imported package to a file.
    """
    with open(platform_file_path, 'a') as file:
        # Write the message about what was just imported
        file.write(f"{message} ")
        
        # Identify the last imported package or component
        if imported_packages:
            full_name = list(imported_packages.keys())[-1]  # Get the last key in the dictionary
            package_name = full_name.split('.')[0]  # Only use the top-level package name for version retrieval
            try:
                if package_name == 'sklearn':
                    package_name = 'scikit-learn'
                pkg_version = version(package_name)
                version_info = f"(version: {pkg_version})"
            except Exception as e:
                version_info = f"(version not found for {package_name}: {e})"
            
            # Write the latest version info to the file and print it
            file.write(f"{version_info}\n\n")
            logging.info(version_info)