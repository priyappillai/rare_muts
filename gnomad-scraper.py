import requests
import urllib.request
import bs4

#url = 'https://gnomad.broadinstitute.org/variant/1-55512261-G-A'
url = 'https://myvariant.info/v1/query?q=CDK7&fields=gnomad_genome&size=1000'
#q = gene name
#fields =  gnomad_genome 
#size = 1000

page = requests.get(url)

if page.status_code == 200:
    print('Success!')
elif page.status_code == 404:
    print('Not Found.')

soup = bs4.BeautifulSoup(page.content, 'lxml')
#table = soup.find(name='table', attrs={'id':'tableID'})

print(soup.prettify())
#a = soup.find_all("div") #, _ngcontent-app-root-c9')
#print(a)

#print(page.content)