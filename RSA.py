import random
import math
import base64
import time

def test_millera_rabina(n, count_round):
	if n == 2 or n == 3 :
		return True
	if (n & 1 == 0) or (n == 1):
		return False
	s = 0
	t = n - 1
	while ((t & 1) != 1):
		t >>= 1
		s += 1
	for i in range(count_round):
		a = random.randint(2, n-2)
		x = pow(a, t, n)
		if (x == 1 or x == n - 1):
			continue
		for j in range(s-1):
			pow(x, 2, n)
			if (x == 1):
				return False
			if (x == n - 1):
				break
		else:
			return False
	return True


def kgen(count_bit):
	
	p = q = 0
	random.seed()
	while p == 0:
		a = random.getrandbits(count_bit-2)
		p = (a<<2)+3
		if (test_millera_rabina(p,int(math.log2(p))) == False):
			p = 0
	while q == 0:
		a = random.getrandbits(count_bit-2)
		q = (a<<2)+3
		if (test_millera_rabina(q, int(math.log2(q))) == False):
			q = 0
	return (p, q, p*q)

def nod_1(a, b):
	if (b == 0):
		return (a, 1, 0)
	(d, kb, kc) = nod_1(b, a%b)
	return (d, kc, kb - (a // b)*kc)



def eratosthenes(b):    
    list_numbers = list(range(b + 1))
    list_numbers[1] = 0    
    for i in list_numbers:
        if i > 1:
            for j in range(i << 1, len(list_numbers), i):
                list_numbers[j] = 0
    return list_numbers

def prime_degree(b, b2):
	list_numbers = eratosthenes(b2)
	j = 0
	k = 0
	mem = 0
	degree_prime = []
	prime_numbers = []
	diff = []
	for i in list_numbers:
		if (i < b):
			if (i != 0):
				degree_prime.append(int(math.log2(b) // math.log2(i)))
				prime_numbers.append(i)
				j += 1
				mem = i
		else:
			diff.append(i - mem)
			mem = i
			k += 1
	return (prime_numbers, degree_prime, diff)

def nod(a, b):
	if (b == 0):
		return (a)
	d = nod(b, a%b)
	return (d)


def pollard_p_1_one(a, n, prime_numbers, degree_prime, time_end):
	x = a
	y = x
	p = 1
	i = j = 0
	d = 1
	i_n = 0
	count = len(prime_numbers)//20 + 1
	for i in range(count):
		if (i != count-1):
			k = 20
		else:
			k = len(prime_numbers) - 20*(count-1)
		for j in range(k):
			x = pow(x, pow(prime_numbers[i*20+j], degree_prime[i*20+j]), n)
			p *= (x-1)
			p %= n
			if (time.time() > time_end):
				return (False, -1)
		d = nod(p, n)
		if (d > 1):
			x = y
			i_n = i
			break
		if (time.time() > time_end):
			return (False, -1)
	if (i == count):
		d = nod(p, n)
		if (d == 1):
			return (False, x)
		i_n = i-1
		x = y
	i = i_n*20
	if (d > 1):
		while (i < len(prime_numbers)):
			for j in range((degree_prime[i])):
				if (time.time() > time_end):
					return (False, -1)
				x = pow(x, prime_numbers[i], n)
				d = nod(x-1, n)
				if (d > 1):
					if (d < n):
						return (True, d)
					else:
						return (False, x)
			i += 1
	return (False, n)

def pollard_p_1(a, n, b1, b2, time_end):
	prime_numbers, degree_prime, diff = prime_degree(b1, b2)
	if (time.time() > time_end):
		return -1
	flag, x = pollard_p_1_one(a, n, prime_numbers, degree_prime, time_end)
	if (x == -1):
		return (False, -1)
	if (flag == True):
		return (True, x)
	b = x
	p = 1
	i = 0
	j = i
	y = x
	mas_b = []
	for i in range(len(diff)):
		mas_b.append(pow(b, diff[i], n))
		if (time.time() > time_end):
			return (False, -1)
	x = pow(x, prime_numbers[len(prime_numbers)-1], n)
	count = len(diff)//20 + 1
	for i in range(count):
		if (time.time() > time_end):
			return (False, -1)
		if (i != count-1):
			k = 20
		else:
			k = len(diff) - 20*(count-1)
		for j in range(k):
			x *= mas_b[i*20 + j]
			x %= n
			p *= (x - 1)
			p %= n
		d = nod(p, n)
		if (d > 1):
			break
		i_n = i 
		y = x
	if (i == count):
		d = nod(p, n)
		if (d == 1):
			return (False, d)
	if (d > 1):
		d = 1
		i = i_n*20+1
		x = y 
		while ((d == 1) and (time.time() < time_end)):
			x *= mas_b[i]
			i += 1
			d = nod(x-1, n)
		if (d < n):
			return (True, d)
	return (False, n)

def p_pollard (n, f, t):
	time_end = time.time() + t
	f.write('Факторизация ро-методом Полларда:\n')
	x = random.randint(2, n-2)
	f.write('    Первый элемент последовательности: ' + str(x))
	i = 0
	y = 1
	state = 2
	while(time.time() < time_end):
		d = nod(x - y, n)
		if ((d > 1) and (d < n) ):
			f.write('    найден множитель: ' + str(d))
			return 0
		if (i == state):
			y = x
			state <<= 1
		x = (x*x + 1) % n
		i+=1
	return -1


def attack_pollard(n, f, t):
	b1 = 100000
	b2 = 1000000
	time_end = time.time() + t
	f.write('Факторизация p-1 методом Полларда:\n')
	flag, d = pollard_p_1(2, n, b1, b2, time_end)
	f.write('    b1 = ' + str(b1) + ' ' + 'b2 =' + str(b2) + '\n')
	f.write('    a = 2 \n')
	if (n == d):
		f.write('    Атака не удалась\n')
		flag, d = pollard_p_1(3, n, b1, b2, time_end)
		f.write('    a = 3 \n')
	if (flag):
		f.write('    Атака прошла успешно, найден множитель: ' + str(d) + '\n')
		return 0
	f.write('    Атака не удалась\n')
	return -1

def chain_fraction(x, y):
	frac = []
	while y:
		a = x // y
		frac.append(a)
		b = x - y*a
		x = y
		y = b
	return frac

def suitable_fraction(frac):
	p_2 = 0
	p_1 = 1
	q_2 = 1
	q_1 = 0
	p_q = []
	for a in frac:
		p = p_1 * a + p_2
		q = q_1 * a + q_2
		p_2 = p_1
		p_1 = p
		q_2 = q_1
		q_1 = q
		p_q.append((p, q))
	return p_q

def attack_vinner(e, n, t, f):
	f.write('Атака Виннера:\n')
	time_end = time.time() + t
	frac = chain_fraction(e, n)
	k_d = suitable_fraction(frac)
	for (k, d) in k_d:
		if (time.time() > time_end):
			return -1
		if (k == 0):
			continue
		phi = (e*d - 1)//k
		dis = math.sqrt(((n-phi+1)*(n-phi+1)-(n*4)))
		p = (n-phi+1 + dis)/2
		q = (n-phi+1 - dis)/2
		if (p*q == n):
			f.write('    Получены следующие множители: ' + str(p) + ' ' + str(q) + '\n')
			return 0
	return -1

def attack_iter_encrypt1(e, n, f):
	random.seed()
	mes = random.randint(2, n)
	cipher = pow(mes, e, n)
	c = pow(cipher, e, n)
	k = 1
	while(c != cipher):
		c1 = c
		c = pow(c, e, n)
		k += 1
		if (k > 80000):
			return (False, n)
	return (True, mes, c1)

def attack_iter_encrypt(e, n, t, f):
	f.write('Итерация процесса шифрования:\n')
	time_end = time.time() + t
	mes = 2
	c1 = 0
	c = pow(mes, e, n)
	d = nod(c - mes, n)
	while(d == 1):
		if (time.time() > time_end):
			return -1
		c = pow(c, e, n)
		d = nod(c - mes, n)
	f.write('    Получен множитель: ' + str(d) + '\n')
	return 0

def div_poly(a, b, time_end):
	dega = len(a)-1
	degb = len(b)-1
	#deg = dega - degb
	#c = a[dega]/b[degb]
	deg = dega - degb
	res = [0]*(deg+1)
	#res[deg] = c
	new_deg = dega
	while (dega >= degb and not(dega == degb == 0 and a[0] == 0)):
		if (time.time() > time_end):
			return -1
		c = a[dega]/b[degb]
		#print(dega)
		deg = dega - degb
		res[deg] = c
		for i in range(dega, -1, -1): 
			if (i - deg >= 0):
				a[i] = a[i] - b[i-deg]*c
			if (a[i] != 0  and new_deg == dega):
				new_deg = i
			if (i == 0 and new_deg == dega):
				new_deg = 0
		dega = new_deg
	q = []
	for i in range(dega+1):
		q.append(a[i])
	#print(res, ' ', q)
	return (q)

def gcd(a,b, time_end):
	dega = len(a) - 1
	degb = len(b) - 1
	while (degb != 0 or b[0] != 0):
		a, b = b, div_poly(a,b, time_end)
		if (time.time() > time_end):
			return -1
		dega = degb
		degb = len(b) - 1
	return a

def factorial(n, p):
	res = 1
	for i in range(1, n+1, 1):
		res *= i
		res %= p
	return res

def nod_inv(a, b):
	if (b == 0):
		return (a, 1, 0)
	(d, kb, kc) = nod_inv(b, a%b)
	return (d, kc, kb - (a // b)*kc)

def inv(a, p):
	d, kc, kb = nod_inv(a,p)
	return kc%p

def binom(a, b, e, p, time_end):
	poly = []
	for i in range(e+1):
		k = factorial(e,p)*inv(factorial(i,p),p)*inv(factorial(e-i,p),p)*pow(b,i,p)*pow(a, n-i,p)
		poly.append(k%p)
		if (time.time() > time_end):
			return -1
	return poly




def attack_lit_deg(e, n, t, f):
	f.write('Атака восстановления открытого текста, использующая маленький модуль шифрования:\n')
	time_end = time.time() + t
	while (time.time() < time_end):
		random.seed()
		m1 = random.randint(2, n)
		a = random.randint(1, 100)
		b = random.randint(1, 100) 
		m2 = a * m1 + b
		c1 = pow(m1, e, n)
		c2 = pow(m1, e, n)
		poly1 = []
		poly1.append(-c1)
		if (time.time() > time_end):
			return -1
		for i in range(e-1):
			poly1.append(0)
			if (time.time() > time_end):
				return -1
		poly1.append(1)
		poly2 = binom(a, b, e, n, time_end)
		if (poly2 == -1):
			return -1
		poly2[0] = (poly2[0] - c2)%n
		nod_poly = gcd(poly1, poly2, time_end)
		if (nod_poly == -1):
			return -1
		deg = len(nod_poly)-1
		if (deg >= 2):
			res = n - inv(len-1, n)*nod_poly[len(nod_poly)-2]
			if (time.time() > time_end):
				return -1
			f.write('    Атака прошла успешно, получен открытый текст: ' + str(res%n) + '\n')
			return 0
	return -1






def rsa_kgen(count_bit):
	
	p = q = 0
	random.seed()
	while p == 0:
		a = random.getrandbits(count_bit-2)
		p = (a<<2)+3
		if (test_millera_rabina(p,int(math.log2(p))) == False):
			p = 0
	while q == 0:
		a = random.getrandbits(count_bit-2)
		q = (a<<2)+3
		if (test_millera_rabina(q, int(math.log2(q))) == False):
			q = 0
	phi = (p-1)*(q-1)
	e = random.randint(int(math.sqrt(p-1)), phi)
	while (nod(e, phi) != 1):
		e = random.randint(int(math.sqrt(p-1)), phi)
	a, b, c  = nod_1(e, phi)
	d = b % phi
	return (e, d, p*q)

def key_p_1_pollard(count_bit, f):
	p = q = 0
	random.seed()
	while p == 0:
		if (count_bit > 7):
			a = random.getrandbits(5)
		else:
			a = random.getrandbits(count_bit-2)
		p = (a<<2)+3
		if (test_millera_rabina(p,int(math.log2(p))) == False):
			p = 0
	while q == 0:
		a = random.getrandbits(count_bit-2)
		q = (a<<2)+3
		if (test_millera_rabina(q, int(math.log2(q))) == False):
			q = 0
	phi = (p-1)*(q-1)
	e = random.randint(int(math.sqrt(p-1)), phi)
	while (nod(e, phi) != 1):
		e = random.randint(int(math.sqrt(p-1)), phi)
	a, b, c  = nod_1(e, phi)
	d = b % phi
	f.write('Открытый ключ: ' + str(e) + ', закрытый ключ: ' + str(d) + ', n: ' + str(p*q) + '\n')
	return 0

def key_viner(count_bit, f):
	p = q = 0
	random.seed()
	while p == 0:
		a = random.getrandbits(count_bit-2)
		p = (a<<2)+3
		if (test_millera_rabina(p,int(math.log2(p))) == False):
			p = 0
	while q == 0:
		a = random.getrandbits(count_bit-2)
		q = (a<<2)+3
		if (test_millera_rabina(q, int(math.log2(q))) == False):
			q = 0
	phi = (p-1)*(q-1)
	d = random.randint(1, int(1/3*math.sqrt(math.sqrt(p*q))))
	while (nod(d, phi) != 1):
		d = random.randint(1, int(1/3*math.sqrt(math.sqrt(p*q))))
	a, b, c  = nod_1(d, phi)
	e = b % phi
	f.write('Открытый ключ: ' + str(e) + ', закрытый ключ: ' + str(d) + ', n: ' + str(p*q) + '\n')
	return 0

def key_lit_mod(count_bit, f):
	p = q = 0
	random.seed()
	while p == 0:
		a = random.getrandbits(count_bit-2)
		p = (a<<2)+3
		if (test_millera_rabina(p,int(math.log2(p))) == False):
			p = 0
	while q == 0:
		a = random.getrandbits(count_bit-2)
		q = (a<<2)+3
		if (test_millera_rabina(q, int(math.log2(q))) == False):
			q = 0
	phi = (p-1)*(q-1)
	e = random.randint(1, int(10))
	while (nod(e, phi) != 1):
		e = random.randint(1, int(10))
	a, b, c  = nod_1(e, phi)
	d = b % phi
	f.write('Открытый ключ: ' + str(e) + ', закрытый ключ: ' + str(d) + ', n: ' + str(p*q) + '\n')
	return 0
	
flag1 = int(input("Выберете режим:\n1 - генерация слабого ключа (могут быть взломаны методом Полига-Хеллмана)\n2 - проверка ключа\n"))
f = open("result.txt", 'w')
if (flag1 == 1):
	nbit = int(input("Укажите размер секретного ключа в битах:\n"))
	number = int(input("1 - p-1 метод Полларда, 2 - атака Виннера, 3 - малый ключ шифрования :\n"))
	if (number == 1):
		key_p_1_pollard(nbit, f)
	if (number == 2):
		key_viner(nbit, f)
	if (number == 3):
		key_lit_mod(nbit, f)
if (flag1 == 2):
	e = int(input("Укажите открытый ключ (e): \n"))
	n = int(input("Укажите открытый ключ (n): \n"))
	f.write('Открытый ключ: ' + str(e) + ' ' + str(n) + '\n')
	t = int(input("Укажите максимальное время работы:\n"))
	flag = int(input("Выберете алгоритм: \n1 - p-1 метод Полларда\n2 - ро-метод Полларда\n3 - маленький модуль шифрования\n4 - атака Винера\n5 - итерация процесса шифрования\n6 - Все вышеперечисленные\n"))
	if (flag == 1):
		if (attack_pollard(n,f,t) == -1):
			f.write('    Разложение не найдено \n')
	if (flag == 2):
		if (p_pollard(n,f,t) == -1):
			f.write('    Разложение не найдено \n')
	if (flag == 3):
		if (attack_lit_deg(e, n, t,f) == -1):
			f.write('    Получить открытый текст не удалось \n')
	if (flag == 4):
		if (attack_vinner(e, n,t, f) == -1):
			f.write('    Разложение не найдено \n')
	if (flag == 5):
		if (attack_iter_encrypt(e, n,t, f) == -1):
			f.write('    Разложение не найдено \n')
	if (flag == 6):
		t = t/5
		if (attack_pollard(n,f,t) == -1):
			f.write('    Разложение не найдено \n')
		if (p_pollard(n,f,t) == -1):
			f.write('    Разложение не найдено \n')
		if (attack_lit_deg(e, n, t,f) == -1):
			f.write('    Получить открытый текст не удалось \n')
		if (attack_vinner(e, n,t, f) == -1):
			f.write('    Разложение не найдено \n')
		if (attack_iter_encrypt(e, n,f,t) == -1):
			f.write('    Разложение не найдено \n')

f.close()









