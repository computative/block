TEMPLATE = app
CONFIG -= console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++17

QMAKE_LIBS += -lstdc++fs
QMAKE_LIBS += -pthread
SOURCES += main.cpp \
    jobmanager.cpp \
    job_info.cpp \
    file_info.cpp \
    buffer.cpp
#QMAKE_CXXFLAGS_RELEASE -= -O
#QMAKE_CXXFLAGS_RELEASE -= -O1
#QMAKE_CXXFLAGS_RELEASE -= -O2
#QMAKE_CXXFLAGS_RELEASE *= -O3

HEADERS += \
    jobmanager.h \
    job_info.h \
    file_info.h \
    buffer.h
