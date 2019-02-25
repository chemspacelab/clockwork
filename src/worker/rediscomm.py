#!/usr/bin/env python
import time
import os
import redis
import sys
import subprocess
import shlex
import gzip
import hashlib

class Taskqueue(object):
	def __init__(self, connectionstring, projectprefix):
		""" Connection string: password@host:6380/1 """
		self._prefix = projectprefix
		password, rest = connectionstring.split('@')
		hostname, rest = rest.split(':')
		port, db = rest.split('/')
		try:
			self._con = redis.Redis(host=hostname, port=int(port), db=int(db), password=password)
		except:
			raise ValueError('Unable to connect to this redis instance.')

	def _start_task(self):
		""" Fetches a new task from redis and returns the associated string."""
		task = self._con.rpoplpush(self._prefix + '_Queue', self._prefix + '_Running')
		if task is None:
			return None

		self._taskid = hashlib.md5(task).hexdigest()
		self._starttime = time.time()
		self._con.hset(self._prefix + '_Started', self._taskid, self._starttime)
		return task
		
	def _store_result(self, taskstring, resultstring, logmessage=None):
		""" Signals task completion to redis and stores the results."""
		# remove server-side backup
		pipeline = self._con.pipeline()
		pipeline.lrem(self._prefix + '_Running', 1, taskstring)
		pipeline.hdel(self._prefix + '_Started', self._taskid)

		# Store results and logs
		pipeline.lpush(self._prefix + '_Results', resultstring)
		if logmessage is not None:
			pipeline.hset(self._prefix + '_Log', taskstring, logmessage)
		pipeline.execute()

		# Update stats
		unixminute = int(time.time() / 60)
		duration = time.time() - self._starttime
		statskey = self._prefix + '_Stats:%d' % unixminute
		pipeline = self._con.pipeline()
		pipeline.rpush(statskey, duration)
		pipeline.expire(statskey, 3600)
		pipeline.execute()
	
	def _get_deadline(self):
		""" Calculates the deadline of job termination. Returns a unix time stamp. When run outside SLURM, returns 0. """
		jobid = os.getenv('SLURM_JOB_ID')
		if jobid is None:
			return 0

		cmd = 'squeue -h -j %s -o "%%L"' % jobid
		try:
			p = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE)
		except:
			return 0
		remaining = p.stdout.decode('utf8').strip()

		# parse time format
		days = 0
		if '-' in remaining:
			days, rest = remaining.split('-')
		else:
			rest = remaining
		try:
			hours, minutes, seconds = rest.split(':')
		except:
			return 0

		return time.time() + ((int(days)*24 + int(hours))*60 + int(minutes))*60 + int(seconds)

	def main_loop(self, callback):
		deadline = self._get_deadline()
		while True:
			if deadline > 0 and deadline - time.time() < 120:
				break
			task_str = self._start_task()
			if task_str is None:
				break
			result_str, log_str = callback(task_str)
			self._store_result(task_str, result_str, log_str)

	def insert(self, taskstring):
		""" Enqueues a task. """
		self._con.lpush(self._prefix + '_Queue', taskstring)

	def get_results(self, purge_after=False):
		""" Fetches and optionally deletes the results."""
		results = self._con.lrange(self._prefix + '_Results')
		results = [gzip.decompress(_).decode('utf8') for _ in results]
		if purge_after:
			self._con.delete(self._prefix + '_Results')
		return results

def do_work(task):
	""" Sample task evaluation. Returns a result as string and an optional log message."""
	time.sleep(2)
	return "result", "Another information"

if __name__ == '__main__':
	tasks = Taskqueue(os.getenv('CHEMSPACELAB_REDIS_CONNECTION'), 'DEBUG')
	for i in range(3):
		tasks.insert('test%d' % i)
	try:
		tasks.main_loop(do_work)
	except:
		# Ugly catch-all handler: no SLURM log file wanted!
		sys.exit(1)
