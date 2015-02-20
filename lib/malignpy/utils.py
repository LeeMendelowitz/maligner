

def is_iterable(obj):
	try:
		i = iter(obj)
	except TypeError:
		return False
	return True