#ifndef PY_MANAGER_HPP
#define PY_MANAGER_HPP

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

inline void initialize_numpy() {
    if (_import_array() < 0) {
        std::cerr << "Failed to initialize NumPy!" << std::endl;
        Py_Finalize();
        return;
    }
}

template <typename T>
int call_python_w_numpy(const T* data, const npy_intp* dims, const std::string& script_name, const std::unordered_map<std::string, std::string>& kwargs = {}) {
    // Initialize the Python interpreter
    Py_Initialize();
    
    // Initialize NumPy
    initialize_numpy(); 

    // Create a NumPy array from the raw data pointer (without deep copy)
    PyObject* pArray = PyArray_SimpleNewFromData(dims ? 1 : 0, dims, NPY_DOUBLE, const_cast<T*>(data));
    if (!pArray) {
        std::cerr << "Error creating NumPy array from data." << std::endl;
        Py_Finalize();
        return -1;
    }

    // Add the script directory to Python's path (if not already included)
    std::string scriptDir = "../scripts"; // Path to the script directory relative to the executable
    PyObject* sysPath = PyImport_ImportModule("sys");
    if (sysPath) {
        PyObject* pPath = PyObject_GetAttrString(sysPath, "path");
        if (pPath) {
            PyList_Append(pPath, PyUnicode_DecodeFSDefault(scriptDir.c_str()));
            Py_XDECREF(pPath);
        }
        Py_XDECREF(sysPath);
    }

    // Import the Python module (e.g., script_name)
    PyObject* pModuleName = PyUnicode_DecodeFSDefault(script_name.c_str());
    PyObject* pModule = PyImport_Import(pModuleName);
    Py_XDECREF(pModuleName);

    if (pModule != nullptr) {
        // Get the function from the Python module
        PyObject* pFunc = PyObject_GetAttrString(pModule, "process_data_from_cpp");

        if (pFunc && PyCallable_Check(pFunc)) {
            // Create kwargs as a Python dictionary if provided
            PyObject* pKwargs = PyDict_New();
            for (const auto& pair : kwargs) {
                PyDict_SetItemString(pKwargs, pair.first.c_str(), PyUnicode_FromString(pair.second.c_str()));
            }

            // Call the Python function, passing the NumPy array and the kwargs (if any)
            PyObject* pArgs = PyTuple_Pack(2, pArray, pKwargs);
            PyObject* pValue = PyObject_CallObject(pFunc, pArgs);
            Py_XDECREF(pArgs);

            if (pValue != nullptr) {
                // If needed, you can process the result from Python here
                Py_XDECREF(pValue);
            } else {
                PyErr_Print(); // Print error if the Python function call failed
            }
        } else {
            std::cerr << "Cannot find function 'process_data_from_cpp' in the Python script." << std::endl;
        }

        Py_XDECREF(pFunc);
        Py_XDECREF(pModule);
    } else {
        PyErr_Print(); // If the Python module failed to load
    }

    // Finalize the Python interpreter
    Py_Finalize();
    return 0;
}


#endif // PY_MANAGER_HPP