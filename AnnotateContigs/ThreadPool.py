# ThreadPool.py - A simple implementation of a threading pool. Gratuitously
# plagiarized this: https://www.metachris.com/2016/04/python-threadpool/
# but with a few modifications for my sanity
#
# Written By: Matt Preston (and Chris Hager)
# Written On: May 31, 2017
# Revised On: Never

try:
    from Queue import Queue
except Exception as e:
    from queue import Queue
from threading import Thread
import logging

class Worker(Thread):
    """Thread executing tasks from a given tasks queue"""
    
    def __init__(self, tasks, name=None):
        Thread.__init__(self, name=name)
        self.tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            try:
                func, args, kargs = self.tasks.get()
                func(*args, **kargs)
            except Exception as e:
                logging.error(e, exc_info=True)
            finally:
                self.tasks.task_done()

class ThreadPool:
    """Pool of threads consuming tasks from a queue"""
    
    def __init__(self, numThreads=0, names=None):
        self.tasks = Queue()
        self.numThreads = 0
        self.AddThreads(numThreads, names)    

    def AddThreads(self, numThreads, names=None):
        """Adds threads to the pool"""
        
        if isinstance(names, list): # Names given
            assert numThreads == len(names), \
                "Unequal pairing of names to number of threads"
            for name in names:
                Worker(self.tasks, name)
        elif names is None:         # No names
            for _ in range(numThreads):
                Worker(self.tasks)
        else:                       # Invalid type for 'names'
            logger = logging.getLogger("%s.%s" % (__name__,
                                                  self.__class__.__name__))
            logger.warning("ThreadPool.__init__(): 'names' is not a list or "
                           "not None, ignoring")
            for _ in range(numThreads):
                Worker(self.tasks)
        self.numThreads += numThreads
            
    def AddTask(self, func, *args, **kargs):
        """Add a task to the queue"""
        
        self.tasks.put((func, args, kargs))

    def Join(self):
        """Wait for completion of all the tasks in the queue"""
        
        self.tasks.join()
