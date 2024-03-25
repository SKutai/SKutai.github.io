import { CommonModule } from '@angular/common';
import { Component, OnDestroy } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Subject } from 'rxjs';
import { takeUntil } from 'rxjs/operators';

@Component({
  selector: 'app-navbar',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './navbar.component.html',
  styleUrl: './navbar.component.scss'
})
export class NavbarComponent implements OnDestroy {
  currentPage: string = 'home';
  private unsubscribe$ = new Subject<void>();

  constructor(private router: Router, private route: ActivatedRoute) {
    this.router.events
      .pipe(takeUntil(this.unsubscribe$))
      .subscribe(() => {
        this.currentPage = this.router.url;
        console.log('Current page is: ', this.currentPage);
      });
  }

  setCurrentPage(page: string) {
    this.currentPage = page;
    console.log('Current page is: ', this.currentPage);
  }

  ngOnDestroy() {
    this.unsubscribe$.next();
    this.unsubscribe$.complete();
  }
}
